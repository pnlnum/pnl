
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

dat = fopen ("complex_bessel_test.dat", mode="w");

function y = besselh1 (nu, x)
    y = besselh (nu, 1, x);
end


function y = besselh2 (nu, x)
    y = besselh (nu, 2, x);
end

Lfunc = { {"besseli", "pnl_complex_bessel_i"},...
    {"besselj", "pnl_complex_bessel_j"},...
    {"besselk", "pnl_complex_bessel_k"},...
    {"bessely", "pnl_complex_bessel_y"},...
    {"besselh1", "pnl_complex_bessel_h1"},...
    {"besselh2", "pnl_complex_bessel_h2"}...
    };

Lfunc_scaled = { {"besseli", "pnl_complex_bessel_i_scaled"},...
    {"besselj", "pnl_complex_bessel_j_scaled"},...
    {"besselk", "pnl_complex_bessel_k_scaled"},...
    {"bessely", "pnl_complex_bessel_y_scaled"}...
};



for j=1:size(Lfunc,2)
    func = Lfunc{j};
    for ii=(1:10)
        arg_r = rand(1) * 20 - 10;
        arg_i = rand(1) * 20 - 10;
        nu = rand(1) * 20 - 10;
        res = feval(func{1}, nu, arg_r + I * arg_i);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f, %.18f},\n", func{2}, func{2}, arg_r, arg_i, nu, real(res), imag(res));
    end
end
for j=1:size(Lfunc_scaled,2)
    func = Lfunc_scaled{j};
    for ii=(1:10)
        arg_r = rand(1) * 20 - 10;
        arg_i = rand(1) * 20 - 10;
        nu = rand(1) * 20 - 10;
        res = feval(func{1}, nu, arg_r + I * arg_i, 1);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f, %.18f},\n", func{2}, func{2}, arg_r, arg_i, nu, real(res), imag(res));
    end
end
fclose(dat)


dat = fopen ("real_besselh_test.dat", mode="w");
for ii=(1:10)
    arg = rand(1) * 20 - 10;
    nu = rand(1) * 20 - 10;
    res = besselh1 (nu, arg);
    fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f},\n", "pnl_bessel_h1", "pnl_bessel_h1", arg, nu, real(res), imag(res));
    res = besselh2 (nu, arg);
    fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f},\n", "pnl_bessel_h2", "pnl_bessel_h2", arg, nu, real(res), imag(res));
end
fclose(dat)

