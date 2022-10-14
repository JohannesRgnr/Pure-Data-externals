/**
 * @file ota~.c
 * @author johannes regnier
 * @brief Simulation of a 4-poles OTA Ladder Filter 
 * @version 0.1
 */


/* based on Miller Puckette's bob~ (Runge-Kutta 4th order)*/
/* copyright 2015 Miller Puckette - BSD license */



#include "m_pd.h"
#include <math.h>
#define DIM 4
#define FLOAT double



typedef struct _params
{
    FLOAT p_input;
    FLOAT p_cutoff;
    FLOAT p_resonance;
    FLOAT p_derivativeswere[DIM];
} t_params;


typedef struct _ota
{
    t_object x_obj;
    t_float x_f;
    t_outlet *x_out1;    /* signal output */

    t_params x_params;
    FLOAT x_state[DIM];
    FLOAT x_sr;
    int x_oversample;
    int x_mode;
    FLOAT p_input;
    FLOAT p_cutoff;
    FLOAT p_resonance;
    FLOAT p_derivativeswere[DIM];

} t_ota;

static void calc_derivatives(FLOAT *dstate, FLOAT *state, t_ota *x)
{
    FLOAT k = ((float)(2*3.14159)) * x->p_cutoff;
   
        dstate[0] = k * tanh(1.1*x->p_input - x->p_resonance*tanh(1.96*state[3]) - state[0]); 
        dstate[1] = k * tanh(1.1*state[0] - state[1]);
        dstate[2] = k * tanh(1.1*state[1] - state[2]);
        dstate[3] = k * tanh(1.1*state[2] - state[3]);
      
}


static void solver_rungekutta(FLOAT *state, FLOAT stepsize, t_ota *x)
{
    FLOAT cumerror = 0;
    int i;
    FLOAT deriv1[DIM], deriv2[DIM], deriv3[DIM], deriv4[DIM], tempstate[DIM];
    FLOAT oldstate[DIM], backstate[DIM];

    calc_derivatives(deriv1, state, x);
    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + 0.5 * stepsize * deriv1[i];
    calc_derivatives(deriv2, tempstate, x);
    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + 0.5 * stepsize * deriv2[i];
    calc_derivatives(deriv3, tempstate, x);
    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + stepsize * deriv3[i];
    calc_derivatives(deriv4, tempstate, x);
    for (i = 0; i < DIM; i++)
        state[i] += (1./6.) * stepsize * 
            (deriv1[i] + 2. * deriv2[i] + 2. * deriv3[i] + deriv4[i]);
}

static t_class *ota_class;


static void ota_oversample(t_ota *x, t_float oversample)
{
    if (oversample <= 1)
        oversample = 1;
    if (oversample > 8)
        oversample = 8;
    x->x_oversample = oversample;
}

static void ota_clear(t_ota *x)
{
    int i;
    for (i = 0; i < DIM; i++)
        x->x_state[i] = x->p_derivativeswere[i] = 0;
}


static void ota_print(t_ota *x)
{
    int i;
    for (i = 0; i < DIM; i++)
        post("state %d: %f", i, x->x_state[i]);
        post("oversample %d", x->x_oversample);
        
}

static void *ota_new( void)
{
    t_ota *x = (t_ota *)pd_new(ota_class);
    x->x_out1 = outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_f = 0;
    ota_clear(x);
    ota_oversample(x, 2);
    return (x);
}

static t_int *ota_perform(t_int *w)
{
    t_ota *x = (t_ota *)(w[1]);
    t_float *in1 = (t_float *)(w[2]);
    t_float *cutoffin = (t_float *)(w[3]);
    t_float *resonancein = (t_float *)(w[4]);
    t_float *out = (t_float *)(w[5]);
    int n = (int)(w[6]), i, j;
    FLOAT stepsize = 1./(x->x_oversample * x->x_sr);
    
    for (i = 0; i < n; i++)
    {
        x->p_input = *in1++;
        x->p_cutoff = *cutoffin++;
        if ((x->p_resonance = *resonancein++) < 0)
            x->p_resonance = 0;
        for (j = 0; j < x->x_oversample; j++)
            solver_rungekutta(x->x_state, stepsize, x);
            *out++ = x->x_state[3]; // Low pass and band pass modes
           }
    return (w+7);
}

static void ota_dsp(t_ota *x, t_signal **sp)
{
    x->x_sr = sp[0]->s_sr;
    dsp_add(ota_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec,
        sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

void ota_tilde_setup(void)
{
    int i;
    ota_class = class_new(gensym("ota~"),
        (t_newmethod)ota_new, 0, sizeof(t_ota), 0, 0);

    class_addmethod(ota_class, (t_method)ota_oversample, gensym("oversample"),
        A_FLOAT, 0);
    class_addmethod(ota_class, (t_method)ota_clear, gensym("clear"), 0);
    class_addmethod(ota_class, (t_method)ota_print, gensym("print"), 0);

    class_addmethod(ota_class, (t_method)ota_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(ota_class, t_ota, x_f);
}
