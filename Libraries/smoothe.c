// Kallyn 2024
// Adapted from lop2~.c from ELSE library, Porres 2017

#include <m_pd.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct _smoothe{
    t_object    x_obj;
    t_int       x_n;
    t_inlet    *x_inlet_freq;
    t_outlet   *x_out;
    t_float     x_nyq;
    double      x_xnm1;
    double      x_ynm1;
    double      x_f;
    double      x_reson;
    double      x_a0;
    double      x_a1;
    double      x_b1;
}t_smoothe;

static t_class *smoothe_class;

static void update_coeffs(t_smoothe *x, double f){
    x->x_f = f;
    double omega = f * PI/x->x_nyq;
    if(omega < 0)
        omega = 0;
    if(omega > 2){
        x->x_a0 = 1;
        x->x_a1 = x->x_b1 = 0;
    }
    else{
        x->x_a0 = x->x_a1 = omega * 0.5;
        x->x_b1 = 1 - omega;
    }
}

static t_int *smoothe_perform(t_int *w){
    t_smoothe *x = (t_smoothe *)(w[1]);
    int nblock = (int)(w[2]);
    t_float *in1 = (t_float *)(w[3]);
    t_float *in2 = (t_float *)(w[4]);
    t_float *out = (t_float *)(w[5]);
    double xnm1 = 0;
    double ynm1 = 0;

    t_int i = 0;

    while (i++ < nblock){
        double xn = *in1++, f = *in2++, yn;
        if(f < 0)
            f = 0;
        else{
            if(f != x->x_f)
                update_coeffs(x, (double)f);
            yn = x->x_a0 * xn + x->x_a1 * xnm1 + x->x_b1 * ynm1;
            *out++ = yn;
            xnm1 = xn;
            ynm1 = yn;
        }
    }

    xnm1 = 0;
    ynm1 = 0;

    while (--i >= 0){
        double xn = *in1++, f = *in2++, yn;
        if(f < 0)
            f = 0;
        else{
            if(f != x->x_f)
                update_coeffs(x, (double)f);
            yn = x->x_a0 * xn + x->x_a1 * xnm1 + x->x_b1 * ynm1;
            *out++ = yn;
            xnm1 = xn;
            ynm1 = yn;
        }
    }

    return(w+6);
}

static void smoothe_dsp(t_smoothe *x, t_signal **sp){
    t_float nyq = sp[0]->s_sr / 2;
    if(nyq != x->x_nyq){
        x->x_nyq = nyq;
        update_coeffs(x, x->x_f);
    }
    dsp_add(smoothe_perform, 5, x, sp[0]->s_n, sp[0]->s_vec,
            sp[1]->s_vec, sp[2]->s_vec);
}

static void smoothe_clear(t_lop2 *x){
    x->x_xnm1 = x->x_ynm1 = 0.;
}

static void *smoothe_new(t_floatarg f){
    t_smoothe *x = (t_smoothe *)pd_new(smoothe_class);
    double freq = f < 0 ? 0 : (double)f;
    x->x_nyq = sys_getsr()/2;
    update_coeffs(x, (double)freq);
    x->x_inlet_freq = inlet_new((t_object *)x, (t_pd *)x, &s_signal, &s_signal);
    pd_float((t_pd *)x->x_inlet_freq, freq);
    x->x_out = outlet_new((t_object *)x, &s_signal);
    return(x);
}

void smoothe_tilde_setup(void){
    smoothe_class = class_new(gensym("smoothe~"), (t_newmethod)smoothe_new, 0,
        sizeof(t_smoothe), CLASS_DEFAULT, A_DEFFLOAT, 0);
    class_addmethod(smoothe_class, (t_method)smoothe_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(smoothe_class, nullfn, gensym("signal"), 0);
    class_addmethod(smoothe_class, (t_method)smoothe_clear, gensym("clear"), 0);
}