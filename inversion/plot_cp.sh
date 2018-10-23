#mkdir plots_inv
#mkdir plots_inv/kerns
#mkdir plots_inv/tests
#mkdir plots_inv/sun
#mkdir plots_inv/phy

cp -R plots/phys/paper/thin/short/inversion-cross-*.pdf plots_inv/kerns
cp -R plots/phys/paper/thin/short/inversion-avg-*.pdf   plots_inv/kerns
cp -R plots/phys-comp/paper/thin/short/*.pdf            plots_inv/phy
cp -R plots/sun/paper/thin/short/*.pdf                  plots_inv/sun
cp -R plots/tests/paper/thin/short/*-modmix.pdf         plots_inv/tests
cp -R plots/tests/paper/thin/short/*-residuals-line.pdf plots_inv/tests
cp -R plots/tests/paper/thin/short/*-solar.pdf          plots_inv/tests
