/* zdsv~ - A zero delay feedback state variable filter */
/* Based on A.Zavalishin TPT*/



#include "m_pd.h"
#include <math.h>
#define FLOAT double


typedef struct _zdsv
{
    t_object x_obj;
    t_float x_f;
    t_outlet *x_out1;   /* LP signal output */
    t_outlet *x_out2;   /* BP signal output */
    t_outlet *x_out3;   /* HP signal output */
    
    FLOAT x_sr;
    FLOAT T; // sampling period
    FLOAT x_lp;
    FLOAT x_bp;
    FLOAT x_hp;
    FLOAT p_input;
    FLOAT p_out;
    FLOAT p_cutoff;
    FLOAT cutoffold;
    FLOAT cutoffincrement;
    FLOAT p_resonance;
    FLOAT resonanceold;
    FLOAT resonanceincrement;
    FLOAT s1;
    FLOAT s2;
    FLOAT PI;

} t_zdsv;



static t_class *zdsv_class;



static void *zdsv_new( void)
{
    t_zdsv *x = (t_zdsv *)pd_new(zdsv_class);
    x->x_out1 = outlet_new(&x->x_obj, gensym("signal"));
    x->x_out2 = outlet_new(&x->x_obj, gensym("signal"));
    x->x_out3 = outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_f = 0;
    x->s1 = 0;
    x->s2 = 0;
    x->PI = 4.0f * atanf(1.0f);
    x->p_cutoff = x->cutoffold = 0.0f;
    x->p_resonance = x->resonanceold = 1.0f;
    return (x);
}

static t_int *zdsv_perform(t_int *w)
{
    t_zdsv *x = (t_zdsv *)(w[1]);
    t_float *in1 = (t_float *)(w[2]);
    t_float *cutoffin = (t_float *)(w[3]);
    t_float *resonancein = (t_float *)(w[4]);
    t_float *out1 = (t_float *)(w[5]);
    t_float *out2 = (t_float *)(w[6]);
    t_float *out3 = (t_float *)(w[7]);
    int n = (int)(w[8]), i, j;
    x->T = 1.0f / x->x_sr; // sampling period
    FLOAT oneoverblocksize = 1.0f/n;
     
    x->p_cutoff = *cutoffin++;
    if(x->p_cutoff != x->cutoffold) // clip cutoff values
    {
       if(x->p_cutoff > x->x_sr*0.48f)
            x->cutoffold = x->p_cutoff = x->x_sr*0.48f;
        else if(x->p_cutoff < 0.0003f)
            x->cutoffold = x->p_cutoff = 0.0003f;
        else
            x->cutoffold = x->p_cutoff;
    }
    x->p_resonance = (1 - 0.01f**resonancein++);
    if(x->p_resonance != x->resonanceold) // clip resonance values
    {
       if(x->p_resonance > 1)
            x->resonanceold = x->p_resonance = 1;
        else if(x->p_resonance < 0.0005f)
            x->resonanceold = x->p_resonance = 0.0005f;
        else
            x->resonanceold = x->p_resonance;
    }

    x->cutoffincrement = (x->p_cutoff - x->cutoffold) * oneoverblocksize;
    x->resonanceincrement = (x->p_resonance - x->resonanceold) * oneoverblocksize;

    for (i = 0; i < n; i++)
    {
        x->p_input = *in1++;
        FLOAT wd = 2*x->PI*x->p_cutoff;
        FLOAT wa = (2.0f * x->x_sr) * tan(wd * x->T * 0.5f);
        FLOAT g = wa * x->T * 0.5f;

        x->x_hp = (x->p_input - 2. * x->p_resonance * x->s1 - g * x->s1 - x->s2) / (1. + 2. * x->p_resonance * g + g * g); 
        x->x_bp = g * x->x_hp + x->s1; 
        x->s1 = g * x->x_hp + x->x_bp; // state update in 1st integrator 
        x->x_lp = g * x->x_bp + x->s2; 
        x->s2 = g * x->x_bp + x->x_lp; // state update in 2nd integrator
        *out1++ = x->x_lp;
        *out2++ = x->x_bp;
        *out3++ = x->x_hp;
        x->p_cutoff += x->cutoffincrement;
        x->p_resonance += x->resonanceincrement;
    }
    return (w+9);
}

static void zdsv_dsp(t_zdsv *x, t_signal **sp)
{
    x->x_sr = sp[0]->s_sr;
    dsp_add(zdsv_perform, 8, x, sp[0]->s_vec, sp[1]->s_vec,
        sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
}

void zdsv_tilde_setup(void)
{
    int i;
    zdsv_class = class_new(gensym("zdsv~"),
        (t_newmethod)zdsv_new, 0, sizeof(t_zdsv), 0, 0);
    class_addmethod(zdsv_class, (t_method)zdsv_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(zdsv_class, t_zdsv, x_f);
}
