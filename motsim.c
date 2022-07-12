/*******************************************
 *                            
 * ルンゲクッタ法による 
 * DCモータシミュレーション       
 * サンプルプログラム
 * 
 * 実行した結果が端末に表示されるので
 * 結果をテキストとして保存しgnuplot等で
 * 可視化する。            
 *                            
********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h> //可変長引数関数を使うために必要
#include <math.h>

//モータの定数
const double Lm = 0.3e-3;//Inductance[H]
const double Rm = 1.2e0;//Resistance[Ohm]1.2e-1
const double Km = 2.3e-2;//Torque constant[Nm/A]
const double Jm = 8.1e-6;//Moment of inertia[kg m^2]
const double Cq = 0.0;//Cofficient of torque (Propeller)
const double Dm = 7.0e-5;   //Cofficient of viscous damping [Nm s]
const double End_time = 1.0;//Time [s]

//モータ構造体
typedef struct 
{
    double i;
    double omega;
    double u;
    double i_;
    double omega_;
    double u_;
} motor_t;

//PWM構造体
typedef struct 
{
    double f;
    double T;
    double duty;
    double high;
    double low;
    double sig;
    double sig_;
    double t_;
} pwm_t;

//フィルター構造体
typedef struct
{
  double tau;
  double h;
  double y;
} filter_t;

//motor構造体の状態変数を保存
void save_motor_state(motor_t *mot)
{
    mot->i_ = mot->i;
    mot->omega_ = mot->omega;
    mot->u_ = mot->u;
}

//モータの状態とフィルターの状態を表示
void print_motor_state(double t, motor_t mot, filter_t filter)
{
    double i=mot.i;
    double omega = mot.omega;
    double u = mot.u_;
    double y = filter.y;
    printf("%11.8f %11.8f %11.8f %11.8f %11.8f\n",t, i, omega, y, u);
}

//PWM構造体を初期化
void pwm_init(pwm_t* pwm, double f, double low, double high, double duty)
{
  pwm->f = f;
  pwm->T = 1/f;
  pwm->low = low;
  pwm->high = high;
  pwm->sig = low;
  pwm->sig_ = low;
  pwm->duty = duty;
}

//PWM構造体のDutyを変更
void pwm_set_duty(pwm_t* pwm, double duty)
{
  pwm->duty = duty;
}

//PWM信号の値を時間を引数にして計算
double pwm_output(pwm_t* pwm, double t)
{
  double out, duty;
  t = fmod(t, pwm->T); 
  duty = pwm->duty * pwm->T;
  if(duty>t) out = pwm->high;
  else out = pwm->low;
  pwm->sig_ = pwm->sig;
  pwm->sig = out;

  return out;
}

//フィルターを初期化
void filter_init(filter_t* filter, double y0, double tau, double h)
{
  filter->y = y0;
  filter->tau = tau;
  filter->h = h;
}

//フィルタの刻み幅セット
void filter_set_h(filter_t* filter, double h)
{
  filter->h = h;
}

//フィルターの更新
void filter_update(filter_t* filter, double u)
{
  double y = filter->y;
  double tau = filter->tau;
  double h = filter->h;
  filter->y = tau*y/(tau+h) + h*u/(tau+h);
}

//Equation of current
//Lm di/dt + Rm i + Km omega = u
//i:current
//t:time
//value[0]:omega
//value[1]:u
double i_dot(double i, double t, double *value)
{
  double omega = value[0];
  double u = value[1];
  double R = value[2];
  return (u - R * i - Km * omega)/Lm;
}

//Equation of motion
//TL = Cq omega^2
//Jm domega/dt + Dm omega + TL + t_f= Km i
//omega:angular velocity
//t:time
//value[0]:i
double omega_dot(double omega, double t, double *value)
{
  double i=value[0];
  double TL=Cq * omega * omega;//プロペラトルク（ドローンシミュレータの名残）
  double omega0 = 1e-8;
  double t_f0 = 0.02; //動摩擦トルク
  double t_f = t_f0*omega/(omega0 + fabs(omega));//摩擦の計算
  return (Km * i - Dm * omega - TL - t_f)/Jm;
}

//Runge Kutta method
//dxdy:derivative
//x:state
//t:time
//h:step size
//n:argument
double rk4(double (*dxdt)(double, double, double*), double x, double t, double h, int n, ...)
{
  va_list args;
  double *value;
  double k1,k2,k3,k4;

  value=(double*)malloc(sizeof(double) * n);
  va_start(args , n);
  for(int i=0;i<n;i++)
  {
    value[i]=va_arg(args, double);
  }
  va_end(args);
  
  k1 = h * dxdt(x, t, value);
  k2 = h * dxdt(x+0.5*h*k1, t+0.5*h, value);
  k3 = h * dxdt(x+0.5*h*k2, t+0.5*h, value);
  k4 = h * dxdt(x+h*k3, t+h, value);

  free(value);
  
  return x+(k1 + k2*2.0 + k3*2.0 + k4)/6;
}

void main(void)
{
  motor_t motor;
  pwm_t pwm;
  filter_t filter;

  double h = 1e-6;//step size  1e-8 on_highz
  double R;
  double pwm_high = 7.4/2;
  double pwm_low = 0.0;

  //PWM周波数
  double pwm_f = 50e3;

  //PWMDuty
  double pwm_duty =0.0;

  //サンプリング周期
  double sample_T = 1e-6;

  double sample =sample_T;
  int N;

  //スイッチングモード on_break
  //このフラグでスイッチングモードを選択
  //0:On-highZ (ohz)
  //1:On-break(ob)
  int on_break=1;
  
  //動作モード変更フラグ duty_vs_omega
  //Duty対角速度を計算するか、モータの過渡応答を計算するか選択
  //0:モータ過渡応答計算
  //1:Duty vs 角速度計算
  int duty_vs_omega = 1;
  
  if (duty_vs_omega==1) 
  {
    N=11;
    pwm_duty=0.0;
  }
  else N=1;

  for (int i=0;i<N;i++)
  {
    double t = 0.0; //time
    motor.omega = 0.0; //angular verocity
    motor.i     = 0.0; //current
    motor.u     = 0.0; //input voltage

    //フィルタ初期化
    filter_init(&filter, 0.0, 0.15, h);
    //PWM初期化
    pwm_init(&pwm, pwm_f, pwm_low, pwm_high, pwm_duty);
    
    //モータ初期状態量表示
    if(duty_vs_omega==0)print_motor_state(t, motor, filter);
    
    //Main loop
    while(t < End_time )    
    {
      //Save state
      save_motor_state(&motor);

      //Control
      motor.u = pwm_output(&pwm, t);

      //Offの処理
      if(on_break==0&&motor.u<(pwm.high+pwm.low)/2)
      {
        //ハイインピーダンス
        R=1e4;
        if( pwm_duty!=0.0)h = 1e-8;
        else h = 1e-7;
      }
      else 
      {
        //オンブレーキ
        R=Rm;
        h = 1e-7;
      }

      //Update(Runge-Kutta method)
      motor.i = rk4(i_dot, motor.i_, t, h, 3, motor.omega_, motor.u, R);
      motor.omega = rk4(omega_dot, motor.omega_, t, h, 1, motor.i_);
      t = t + h;

      //フィルター処理
      filter_set_h(&filter, h);
      filter_update(&filter, motor.omega);

      //Output(sampling)
      if (duty_vs_omega==0 && t>=sample)
      {
        print_motor_state(t, motor, filter);
        sample = sample + sample_T;
      }
      
    }
    //Duty 対　角速度の計算結果を表示
    if (duty_vs_omega==1)
    {
      printf("%11.2f %11.2f %11.2f\n", pwm_duty*100, filter.y, motor.omega);
      fprintf(stderr, "%11.2f %11.2f %11.2f\n", pwm_duty*100, filter.y, motor.omega);
    }
    pwm_duty = pwm_duty + 0.1;
  }
}
