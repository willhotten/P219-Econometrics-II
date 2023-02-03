%%% OLS function

function result = ols(y,X)
    result.beta = (X'*X)\(X'*y);
end