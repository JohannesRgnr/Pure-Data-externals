/* ring64~ - A zero delay feedback 64 bands resonator */
/* Bandpass filter based on A.Zavalishin "The Art of VA Filter Design 1.1.1" and on code from Robin Schmidt*/

/*  Johannes Regnier 04/2018 jregnier@ucsd.edu */

#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#define FLOAT double
#define BANDS 64

typedef struct _ring64
{
    t_object x_obj;
    t_float x_f;
    t_outlet *x_out;   // main signal output

    FLOAT x_sr; // sampling frequency
    FLOAT T; // sampling period
    FLOAT PI;
    
    /* input parameters */
    FLOAT p_input;
    FLOAT p_cutoff;
    FLOAT cutoffold;
    FLOAT cutoffincrement;
    FLOAT p_resonance;
    FLOAT resonanceold;
    FLOAT resonanceincrement;
    FLOAT p_gainband[BANDS];
    FLOAT gainbandold[BANDS];
    FLOAT gainbandincrement[BANDS];
    FLOAT p_brightness;
    FLOAT brightnessold;
    FLOAT brightnessincrement;

    t_int numberbands; // number of active bands
    t_int softclip;

    FLOAT x_lp[BANDS];
    FLOAT x_bp[BANDS];
    FLOAT x_hp[BANDS];
    FLOAT s1[BANDS];
    FLOAT s2[BANDS];

    FLOAT freqmult[BANDS];
    FLOAT freqband[BANDS];    
    FLOAT gainband[BANDS];
    FLOAT singleout[BANDS];
    FLOAT sumout;
    FLOAT gain; // main gain
} t_ring64;



static t_class *ring64_class;



static void *ring64_new( void)
{
    t_ring64 *x = (t_ring64 *)pd_new(ring64_class);
    x->x_out = outlet_new(&x->x_obj, gensym("signal"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);    
    x->x_f = 0;
    x->PI = 4.0f * atanf(1.0f);

    /* default values */
    int m;
    for (m = 0; m < BANDS; ++m) // init all 64 filters
    {
        x->s1[m] = 0;
        x->s2[m] = 0;
        x->freqmult[m] = 1;
        x->gainband[m] = 1;
    }
    x->numberbands = 16;
    x->p_cutoff = x->cutoffold = 0.0f;
    x->p_resonance = x->resonanceold = 1.0f;
    x->p_brightness = x->brightnessold = 0;
    x->softclip = 0;
    x->gain = 0.9;
    return (x);
}


void ring64_freqs(t_ring64 *x, t_symbol *selector, int argcount, t_atom *argvec)
{
    int i;
    for (i = 0; i < argcount; i++)
    {
    	if (argvec[i].a_type == A_FLOAT)
      {
        x->freqmult[i]= argvec[i].a_w.w_float;
      }
    	else if (argvec[i].a_type == A_SYMBOL)
	    error("Wrong argument type: %s", argvec[i].a_w.w_symbol->s_name);
    }
}

void ring64_gains(t_ring64 *x, t_symbol *selector, int argcount, t_atom *argvec)
{
    int i;
    for (i = 0; i < argcount; i++)
    {
    	if (argvec[i].a_type == A_FLOAT)
      {
        if(argvec[i].a_w.w_float > 16.0f )
             x->p_gainband[i] = argvec[i].a_w.w_float  = 16.0f;
        else if(argvec[i].a_w.w_float < 0.0f)
             x->p_gainband[i] = argvec[i].a_w.w_float  = 0.0f;
        else x->p_gainband[i]= argvec[i].a_w.w_float;
      }
    	else if (argvec[i].a_type == A_SYMBOL)
	    error("Wrong argument type: %s", argvec[i].a_w.w_symbol->s_name);
    }
}

void ring64_bands(t_ring64 *x, t_float bands)
{
  if(bands<1)
      x->numberbands = 1;
  else if (bands>64)
      x->numberbands = 64;
  else x->numberbands = bands;
}

void ring64_gain(t_ring64 *x, t_float gain)
{
    if (gain<0)
        x->gain = 0;
    else if (gain>2)
            x->gain = 2;
    else x->gain = gain;
}


void ring64_softclip(t_ring64 *x, t_float softclip)
{
  x->softclip = softclip;
}


static void ring64_print(t_ring64 *x)
{
    post("%d bands ", x->numberbands);
    if (x->softclip == 1)
    {
        post("soft clip ON");
    }
    else
        post("soft clip OFF");
}


static t_int *ring64_perform(t_int *w)
{
    t_ring64 *x = (t_ring64 *)(w[1]);
    t_float *in1 = (t_float *)(w[2]);
    t_float *cutoffin = (t_float *)(w[3]);
    t_float *resonancein = (t_float *)(w[4]);
    t_float *p_brightnessin = (t_float *)(w[5]);
    t_float *out = (t_float *)(w[6]);
    int n = (int)(w[7]), i, j;
    x->T = 1.0f / x->x_sr; // sampling period
    FLOAT oneoverblocksize = 1.0f/n;
    FLOAT oneovernumberbands = 1.0f/x->numberbands;


    // if (x->gain<0)
    //     x->gain = 0;
    //     else if (x->gain>2)
    //         x->gain = 2;
           
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
    
    x->p_resonance = 1-exp((-1000.0f/ x->x_sr) / (6.91**resonancein++)); // exponential decay time (empirical)
    if(x->p_resonance != x->resonanceold) // clip resonance values
    {
       if(x->p_resonance > 1)
            x->resonanceold = x->p_resonance = 1;
        else if(x->p_resonance < 0.00002f)
            x->resonanceold = x->p_resonance = 0.00002f;
        else
            x->resonanceold = x->p_resonance;
    }

    x->p_brightness = *p_brightnessin++;
    if (x->p_brightness != x->brightnessold)
    {
        if(x->p_brightness<-1)
            x->brightnessold = x->p_brightness = -1;
        else if (x->p_brightness>1)
            x->brightnessold = x->p_brightness = 1;
        x->brightnessold = x->p_brightness;
    }

 
    x->cutoffincrement = (x->p_cutoff - x->cutoffold) * oneoverblocksize;
    x->resonanceincrement = (x->p_resonance - x->resonanceold) * oneoverblocksize;
    x->brightnessincrement = (x->p_brightness - x->brightnessold) * oneoverblocksize;
    FLOAT pivot = 4; // completely empirical.... could maybe be user defined.. 


    for (i = 0; i < n; i++)
    {

        x->p_input = *in1++;
        
        int k;
        for (k = 0; k < BANDS; ++k)
        {
            if (x->p_cutoff * x->freqmult[k] > 0.48 * x->x_sr) // band limiting
            {
                x->gainband[k] = 0; // mute band
                x->s1[k] = 0; // reset filter states
                x->s2[k] = 0;
            }
            else
            {            
            x->gainband[k] = x->p_gainband[k]*(x->p_brightness*((k+1-pivot)/pivot)+ 1);
            }

            if (x->gainband[k] < 0)
                x->gainband[k] = 0;
        }

        int m;
        for (m = 0; m < x->numberbands; ++m)
        {
            FLOAT wd = 2*x->PI*x->p_cutoff*x->freqmult[m];
            FLOAT wa = (2.0f * x->x_sr) * tan(wd * x->T * 0.5f);
            FLOAT g = wa * x->T * 0.5f;
            x->x_hp[m] = (x->p_input - 2. * x->p_resonance * x->s1[m] - g * x->s1[m] - x->s2[m]) / (1. + 2. * x->p_resonance * g + g * g); 
            x->x_bp[m] = g * x->x_hp[m] + x->s1[m]; 
            x->s1[m] = g * x->x_hp[m] + x->x_bp[m]; // state update in 1st integrator 
            x->x_lp[m] = g * x->x_bp[m] + x->s2[m]; 
            x->s2[m] = g * x->x_bp[m] + x->x_lp[m]; // state update in 2nd integrator
            x->singleout[m] = x->x_bp[m]*x->gainband[m];
        }
        int l;
        x->sumout = 0;
        for (l = 0; l < x->numberbands; ++l)
          {
            x->sumout = x->sumout + (x->gain*x->singleout[l] * oneovernumberbands);
          }
        
        //  soft-clipping if desired.
        if(x->softclip==1)
                {
                    // Limit signal from -1 to 1
                    if (x->sumout > 1.0f)
                        x->sumout = 1.0f;
                    if (x->sumout < -1.0f)
                        x->sumout = -1.0f;
                    *out++ = (1.5f * x->sumout - 0.5f * x->sumout * x->sumout * x->sumout);
                }
            else  *out++ = x->sumout;

        x->p_cutoff += x->cutoffincrement;
        x->p_resonance += x->resonanceincrement;
        x->p_brightness += x->brightnessincrement;

    }
    return (w+8);
}

static void ring64_dsp(t_ring64 *x, t_signal **sp)
{
    x->x_sr = sp[0]->s_sr;
    dsp_add(ring64_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec,
        sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n);
}

void ring64_tilde_setup(void)
{
    int i;
    ring64_class = class_new(gensym("ring64~"),
        (t_newmethod)ring64_new, 0, sizeof(t_ring64), 0, 0);
    class_addmethod(ring64_class, (t_method)ring64_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(ring64_class, (t_method)ring64_freqs, gensym("freqs"), A_GIMME, 0);
    class_addmethod(ring64_class, (t_method)ring64_gains, gensym("gains"), A_GIMME, 0);
    class_addmethod(ring64_class, (t_method)ring64_bands, gensym("bands"), A_FLOAT, 0);
    class_addmethod(ring64_class, (t_method)ring64_print, gensym("print"), 0);
    class_addmethod(ring64_class, (t_method)ring64_softclip, gensym("softclip"), A_FLOAT, 0);
    class_addmethod(ring64_class, (t_method)ring64_gain, gensym("gain"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(ring64_class, t_ring64, x_f);
}
