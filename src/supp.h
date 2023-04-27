#ifndef SUPP_H
#define SUPP_H

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"
#include "random_16807.h"

int PDF2Prob(hist_p_t supp, pdf_p_t pdf, prob_p_t *prob_p);

double SampleRadius(prob_p_t prob, int type_id);

int PDF2Hist(prob_p_t prob, hist_p_t hist);

int Hist2Supp(hist_p_t hist, tbf_p_t *supp_p);

#endif
