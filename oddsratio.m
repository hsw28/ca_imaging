function [h p stats] = oddsratio(vector1, vector2, cutoff1, cutoff2)

%put in the two vectors. any value ABOVE the cutoff will be considered a positive
%returns the p value and odds ratio for the positives in the two vectors

want1 = find(vector1 >= cutoff1);
want2 = find(vector2 >= cutoff2);

holdvalues1 = zeros(length(vector1),1);
holdvalues2 = zeros(length(vector2),1);

holdvalues1(want1) = 1;
holdvalues2(want2) = 1;

CT = crosstab(holdvalues1,holdvalues2);

[h p stats] = fishertest(CT)
