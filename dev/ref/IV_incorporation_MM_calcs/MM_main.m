CFpoints_in = [0; 0.3; 0.6; 0.9]

dat_m6A_T7pol = [ 0.0000484032975; 0.05356719; 0.261309278040897; 0.843465660793357]
dat_m6A_RDpol = [ 0.000160407476252318; 0.0523586599524652; 0.171413407075576; 0.678020892780635  ]

dat_2p0mG_T7    = [0 ; 0.0027918255 ; 0.0049324326; 0.0119157514]
dat_2p0mG_RDpol = [ 0; 0.0024297299; 0.004943885; 0.0349503456 ]

dat_2p0mU_T7 = [0; 0; 0.0019584739; 0.0091498189 ]

% functions defined in adjoinging files

fun_m6A_T7 = @(Pr) MF_model_se(Pr, CFpoints_in, dat_m6A_T7pol);
fun_m6A_RD = @(Pr) MF_model_se(Pr, CFpoints_in, dat_m6A_RDpol);

fun_2p0mG_T7 = @(Pr) MF_model_se(Pr, CFpoints_in, dat_2p0mG_T7);
fun_2p0mG_RD = @(Pr) MF_model_se(Pr, CFpoints_in, dat_2p0mG_RDpol);

fun_2p0mU_T7 = @(Pr) MF_model_se(Pr, CFpoints_in, dat_2p0mU_T7) ;


Pr0 = 1;
% now perform the optimization:

bestPr_m6A_T7 = fminsearch(fun_m6A_T7 , Pr0);
bestPr_m6A_RD = fminsearch(fun_m6A_RD , Pr0);

bestPr_2p0mG_T7 = fminsearch( fun_2p0mG_T7, Pr0);
bestPr_2p0mG_RD = fminsearch( fun_2p0mG_RD, Pr0);

bestPr_2p0mU_T7 = fminsearch( fun_2p0mU_T7, Pr0);

% now plot:
xsample = linspace(0,1, 100)';

% ---- m6A -----
yT_m6A_T7 = MF(bestPr_m6A_T7, xsample);
yT_m6A_RD = MF(bestPr_m6A_RD, xsample);


figure(1)
hold on
plot( CFpoints_in, dat_m6A_T7pol, "o", "color", "blue", 'MarkerSize', 10, 'LineWidth', 4)
plot( xsample, yT_m6A_T7, "--", "color", "blue", 'LineWidth', 4 )

plot( CFpoints_in, dat_m6A_RDpol, "+", "color", "red", 'MarkerSize', 10, 'LineWidth', 4 )
plot( xsample, yT_m6A_RD, "-.", "color", "red", 'LineWidth', 4 )

xlabel ("Input fractioin")
ylabel ("measured output")
title ("m6A")
% ---- 2p0mG -----
yT_2p0mG_T7 = MF(bestPr_2p0mG_T7, xsample);
yT_2p0mG_RD = MF(bestPr_2p0mG_RD, xsample);

figure(2)
hold on
plot( CFpoints_in, dat_2p0mG_T7, "o", "color", "blue")
plot( xsample, yT_2p0mG_T7, "--", "color", "blue" )

plot( CFpoints_in, dat_2p0mG_RDpol, "+", "color", "red")
plot( xsample, yT_2p0mG_RD, "-.", "color", "red" )


