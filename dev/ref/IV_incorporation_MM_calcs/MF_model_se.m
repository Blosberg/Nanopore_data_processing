function result = MF_model_se (Pr, CFpoints, Inc_exp )


 model_prediction = MF (Pr, CFpoints );

result = sum(( model_prediction - Inc_exp).^2);

end