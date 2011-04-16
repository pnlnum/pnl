
// 
//   Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU Lesser General Public License as 
//   published by  the Free Software Foundation; either version 3 of the
//   License, or (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Lesser General Public License for more details.
//  
//   You should have received a copy of the GNU Lesser General Public
//   License  along with this program.  If not, see
//   <http://www.gnu.org/licenses/>. 
//  

dat = mopen ("real_bessel_test.dat", mode="w");
Lfunc = list( ["besseli", "pnl_bessel_i"],...
    ["besselj", "pnl_bessel_j"],...
    ["besselk", "pnl_bessel_k"],...
    ["bessely", "pnl_bessel_y"]...
    );

Lfunc_scaled = list( ["besseli", "pnl_bessel_i_scaled"],...
    ["besselj", "pnl_bessel_j_scaled"],...
    ["besselk", "pnl_bessel_k_scaled"],...
    ["bessely", "pnl_bessel_y_scaled"]...
    );

for func=Lfunc
    for i=(1:10)
        arg = grand(1,1,'unf', -10, 10);
        nu = grand(1,1,'unf', -10, 10);
        res = evstr(func(1) + "(" + string(nu) + "," + string(arg)  + ")"); 
        mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func(2), func(2), arg, nu, res)
    end
end
for func=Lfunc_scaled
    for i=(1:10)
        arg = grand(1,1,'unf', -10, 10);
        nu = grand(1,1,'unf', -10, 10);
        res = evstr(func(1) + "(" + string(nu) + "," + string(arg)  + ",1)"); 
        mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f },\n", func(2), func(2), arg, nu, res)
    end
end

mclose(dat)
