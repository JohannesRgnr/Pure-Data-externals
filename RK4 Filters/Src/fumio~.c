/**
 * @file fumio~.c
 * @author johannes regnier
 * @brief fumio~ - Simulation of a 2-poles Korg MS20 Filter
 * @version 0.1
 * 
 */

/* based on Miller Puckette's bob~ (Runge-Kutta 4th order)*/
/* copyright 2015 Miller Puckette - BSD license */


#include "m_pd.h"
#include <math.h>
#define DIM 2
#define FLOAT double



typedef struct _fumio
{
    t_object x_obj;
    t_float x_f;
    t_outlet *x_out1;    /* signal output */

    // t_params x_params;
    FLOAT x_state[DIM];
    FLOAT x_sr;
    int x_oversample;
    int x_mode;
    FLOAT p_input;
    FLOAT p_cutoff;
    FLOAT p_resonance;
    FLOAT p_derivativeswere[DIM];

} t_fumio;

static void calc_derivatives(FLOAT *dstate, FLOAT *state, t_fumio *x)
{
    FLOAT k = ((float)(2*3.14159)) * x->p_cutoff;


    if (x->x_mode == 1) // low pass
    {
        dstate[0] = k * (x->p_input - state[0] - tanh(x->p_resonance * state[1])); 
        dstate[1] = k * (state[0] - state[1] + tanh(x->p_resonance * state[1]));
    }
    else if (x->x_mode == 3) // high pass
    {
        dstate[0] = k * (state[0] - tanh(x->p_resonance * state[1])); 
        dstate[1] = k * (-x->p_input - state[1]);
    }
    else if   (x->x_mode == 2) // band pass
    {
        dstate[0] = k * ( -x->p_input - state[0] - tanh(x->p_resonance * state[1])); 
        dstate[1] = k * ( x->p_input + state[0] - state[1]+ tanh(x->p_resonance * state[1]));
    }

  
}


static void solver_rungekutta(FLOAT *state, FLOAT stepsize, t_fumio *x)
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

static t_class *fumio_class;


static void fumio_oversample(t_fumio *x, t_float oversample)
{
    if (oversample <= 1)
        oversample = 1;
    if (oversample > 8)
        oversample = 8;
    x->x_oversample = oversample;
}

static void fumio_clear(t_fumio *x)
{
    int i;
    for (i = 0; i < DIM; i++)
        x->x_state[i] = x->p_derivativeswere[i] = 0;
}

static void fumio_mode(t_fumio *x, t_float mode)
{
    int i;
    if (mode >= 1 && mode <=3)
    {
        for (i = 0; i < DIM; i++)
        x->x_state[i] = x->p_derivativeswere[i] = 0;
        x->x_mode = mode; 
    }
    else 
    {
        mode = x->x_mode;
    }      
}

static void fumio_print(t_fumio *x)
{
    int i;
    if (x->x_mode == 1) 
        post("mode: %s", "low pass");
    else if (x->x_mode == 3) 
        post("mode: %s", "high pass");
    else if   (x->x_mode == 2) 
        post("mode: %s", "band pass");
    for (i = 0; i < DIM; i++)
        post("state %d: %f", i, x->x_state[i]);
        post("oversample %d", x->x_oversample);
        
}

static void *fumio_new( void)
{
    t_fumio *x = (t_fumio *)pd_new(fumio_class);
    x->x_out1 = outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_f = 0;
    fumio_clear(x);
    fumio_oversample(x, 2);
    fumio_mode(x, 1);
    return (x);
}

static t_int *fumio_perform(t_int *w)
{
    t_fumio *x = (t_fumio *)(w[1]);
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
        if (x->x_mode == 1 || x->x_mode == 2) 
            *out++ = x->x_state[1]; // Low pass and band pass modes
        else if (x->x_mode == 3) 
            *out++ = x->x_state[1]+x->p_input; // High pass mode
    }
    return (w+7);
}

static void fumio_dsp(t_fumio *x, t_signal **sp)
{
    x->x_sr = sp[0]->s_sr;
    dsp_add(fumio_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec,
        sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

void fumio_tilde_setup(void)
{
    int i;
    fumio_class = class_new(gensym("fumio~"),
        (t_newmethod)fumio_new, 0, sizeof(t_fumio), 0, 0);

    class_addmethod(fumio_class, (t_method)fumio_oversample, gensym("oversample"),
        A_FLOAT, 0);
    class_addmethod(fumio_class, (t_method)fumio_mode, gensym("mode"), A_FLOAT, 0);
    class_addmethod(fumio_class, (t_method)fumio_clear, gensym("clear"), 0);
    class_addmethod(fumio_class, (t_method)fumio_print, gensym("print"), 0);

    class_addmethod(fumio_class, (t_method)fumio_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(fumio_class, t_fumio, x_f);
}
