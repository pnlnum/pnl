
% 
%   Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
%  
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License as 
%   published by  the Free Software Foundation; either version 3 of the
%   License, or (at your option) any later version.
%  
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%  
%   You should have received a copy of the GNU Lesser General Public
%   License  along with this program.  If not, see
%   <http:%www.gnu.org/licenses/>. 
%  

dat = fopen ("gamma_test.dat", mode="w");
Lfunc = { {"gamma", "pnl_sf_gamma"},...
    {"erf", "pnl_sf_erf"}, ...
    {"erfc", "pnl_sf_erfc"}, ...
    {"gammaln", "pnl_sf_log_gamma"}...
    };

for j=1:size(Lfunc,2)
    func = Lfunc{j};
    for i=(1:10)
        arg = rand(1) * 10;
        res = feval(func{1}, arg);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f },\n", func{2}, func{2}, arg, res)
    end
end

fclose(dat)

% Incomplete Gamma function tests
dat = fopen ("gammainc_test.dat", mode="w");

function y = gammainc2 (x, a)
    y = gammainc (x, a, "upper") * gamma(a);
end

function y = gammaincP (x, a)
    y = (1 - gammainc (x, a, "upper"));
end

function y = gammaincQ (x, a)
    y = gammainc (x, a, "upper");
end

function y = choose (n, p)
    y = gamma (n+1) / (gamma (p+1) * gamma (n-p+1));
end

Lfunc = { {"gammainc2", "pnl_sf_gamma_inc"},...
    {"gammaincP", "pnl_sf_gamma_inc_P"}, ...
    {"gammaincQ", "pnl_sf_gamma_inc_Q"}, ...
    };

for j=1:size(Lfunc,2)
    func = Lfunc{j};
    for i=(1:10)
        arg = rand(1) * 10;
        alpha = rand(1) * 10;
        res = feval(func{1}, arg, alpha);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func{2}, func{2}, arg, alpha, res)
    end
end

func = {"choose", "wrap_pnl_sf_choose"};
for i=(1:10)
    n = floor(rand(1) * 30);
    k = floor(rand(1) * n);
    res = feval (func{1}, n, k);
    fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func{2}, func{2}, k, n , res)
end


fclose(dat)

