Bayesian Project References 

- K&P Section 2.4 (censoring mechanism), Sections 7.1–7.3 (competing risks formulation/cause-specific/ subdistribution functions), K&P Section 8.4 (cont. frailty models)
- Hougaard's Analysis of Multivariate Survival Data (2000) comprehensive reference (frailty models). 
- Lancaster (1979, Econometrica) & Vaupel, Manton & Stallard (1979) (connection between frailty and informative censoring) unobserved heterogeneity and selection bias from censoring are fundamentally linked.
- sensitivity analysis:  Scharfstein, Rotnitzky & Robins (1999, JASA)  nonignorable censoring, (formalizes the selection model approach), Rotnitzky, Robins & Scharfstein (1998) extends this to semiparametric settings
- Bayesian treatment: Daniels & Hogan's Missing Data in Longitudinal Studies (2008), pattern-mixture/selection model (directly maps to problem), (Chapter 9 on sensitivity analysis for dropout)
- PG augmentation (Polson, Scott & Windle (2013, JASA) is your core reference. Multinomial extension relevant for competing risks, Linderman, Johnson & Adams (2015) extend the PG trick to multinomial/softmax models. For the practical implementation, the BayesLogit and pgdraw R/Python packages are essential.
**For posterior predictive checks as a censoring diagnostic, Gelman, Meng & Stern (1996, Statistica Sinica) on PPCs generally, and then your specific strategy of comparing against Kaplan-Meier stratified by discharge-time quartile is something you'll need to develop yourself — I'm not aware of a reference that does exactly this, which is actually a strength for the innovation component of your paper.**

