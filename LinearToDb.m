function [ X_abs_Db ] = LinearToDb( X_abs )
%UNTITLED Summary of this function goes here
%   If X(i) > 1 , than X_db(i) = 20*log(X), else X_db(i) = 0
X_abs_Db = zeros(size(X_abs));
X_abs_Db(X_abs > 1) = 20*log10(X_abs(X_abs > 1));

end

