function d = cohens_d(x,y,DIM,NO_SIGN)
% 
% Function will compute cohen's d effect size.
%
% INPUT
%	x = matrix or vector
%	y = matrix or vector
%	DIM = specify the dimension which samples are along
%   NO_SIGN = set to 1 to output effect size ignoring sign, 0 to output with the sign retained
%
% Example usage
%
% x = [1:100]'
% y = [101:200]'
% d = cohens_d(x, y, 1)
%
% written by mvlombardo - 30.08.2015
%

% n-1 for x and y
lx = size(x,DIM)-1;
ly = size(y,DIM)-1;

% mean difference (numerator)
if NO_SIGN
    md = abs(nanmean(x,DIM) - nanmean(y,DIM));
else
    md = nanmean(x,DIM) - nanmean(y,DIM);
end

% pooled variance (denominator)
csd = (lx .* nanvar(x,0,DIM)) + (ly .* nanvar(y,0,DIM));
csd = sqrt(csd./(lx + ly));

% Cohen's d
d = md./csd;

end % function d = cohens_d(x,y,DIM)