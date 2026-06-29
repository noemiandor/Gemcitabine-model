library(xtable)
p=c("$\\iota$","u", "w", "v", "$\\Gamma_{\\alpha}$", "$\\Gamma_a$",  "$\\Gamma_k$", "$\\gamma_a$","$\\lambda$", "$\\mu$", "$\\theta$", "$\\nu$", "$\\eta$", "$\\xi$")
d=c("Maximum number of WGD events before cell enters quiescent state","Drug concentration at which proliferation is half maximum", "Drug concentration at which apoptosis is half maximal", "Drug concentration at which WGD rate is half maximal", "maximum proliferation rate", "maximum apoptosis rate", "maximum WGD rate", "background apoptosis rate", "ploidy conferred decline in drug-induced apoptosis", "ploidy conferred decline/increase in WGD rate","Maximum rate at which Gemcitabine is transported into cells","Gemcitabine concentration at which cell-transporter are at half maximal saturation", "Rate at which Gemcitabine is converted into dFdCTP inside a cell","Rate at which dFdCTP is incorporated into DNA")

tab=as.data.frame(d)
rownames(tab)=p
# colnames(tab)=c("Description")
print(xtable(tab),  include.rownames=T, include.colnames=F, floating=F, sanitize.rownames.function = identity,file = "~/Downloads/Gemcitabinemodelparams.tex")

