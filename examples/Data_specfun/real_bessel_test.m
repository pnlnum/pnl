
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

dat = fopen ("real_bessel_test.dat", mode="w");
Lfunc = { {"besseli", "pnl_bessel_i"},...
    {"besselj", "pnl_bessel_j"},...
    {"besselk", "pnl_bessel_k"},...
    {"bessely", "pnl_bessel_y"}...
};

Lfunc_scaled = { {"besseli", "pnl_bessel_i_scaled"},...
    {"besselj", "pnl_bessel_j_scaled"},...
    {"besselk", "pnl_bessel_k_scaled"},...
    {"bessely", "pnl_bessel_y_scaled"}...
};

for j=1:size(Lfunc,2)
    func = Lfunc{j};
    for i=(1:10)
        arg = rand(1) * 20 - 10;
        nu = rand(1) * 20 - 10;
        res = feval(func{1}, nu, arg);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func{2}, func{2}, arg, nu, res)
    end
end

for j=1:size(Lfunc_scaled,2)
    func = Lfunc_scaled{j};
    for i=(1:10)
        arg = rand(1) * 20 - 10;
        nu = rand(1) * 20 - 10;
        res = feval(func{1}, nu, arg, 1);
        fprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func{2}, func{2}, arg, nu, res)
    end
end
fclose(dat)
