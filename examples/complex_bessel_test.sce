
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

dat = mopen ("complex_bessel_test.dat", mode="w");

function y = besselh1 (nu, x)
    y = besselh (nu, 1, x)
endfunction


function y = besselh2 (nu, x)
    y = besselh (nu, 2, x)
endfunction

Lfunc = list( ["besseli", "pnl_complex_bessel_i"],...
    ["besselj", "pnl_complex_bessel_j"],...
    ["besselk", "pnl_complex_bessel_k"],...
    ["bessely", "pnl_complex_bessel_y"],...
    ["besselh1", "pnl_complex_bessel_h1"],...
    ["besselh2", "pnl_complex_bessel_h2"]...
    );

Lfunc_scaled = list( ["besseli", "pnl_complex_bessel_i_scaled"],...
    ["besselj", "pnl_complex_bessel_j_scaled"],...
    ["besselk", "pnl_complex_bessel_k_scaled"],...
    ["bessely", "pnl_complex_bessel_y_scaled"]...
    );



for func=Lfunc
    for i=(1:10)
        arg_r = grand(1,1,'unf', -10, 10);
        arg_i = grand(1,1,'unf', -10, 10);
        nu = grand(1,1,'unf', -10, 10);
        res = evstr(func(1) + "(" + string(nu) + "," + string(arg_r) + "+ %i *" + string(arg_i)  + ")"); 
        mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f, %.18f},\n", func(2), func(2), arg_r, arg_i, nu, real(res), imag(res));
    end
end
for func=Lfunc_scaled
    for i=(1:10)
        arg_r = grand(1,1,'unf', -10, 10);
        arg_i = grand(1,1,'unf', -10, 10);
        nu = grand(1,1,'unf', -10, 10);
        res = evstr(func(1) + "(" + string(nu) + "," + string(arg_r) + "+ %i *" + string(arg_i)  + ",1)"); 
        mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f, %.18f},\n", func(2), func(2), arg_r, arg_i, nu, real(res), imag(res));
    end
end
mclose(dat)


dat = mopen ("real_besselh_test.dat", mode="w");
for i=(1:10)
    arg = grand(1,1,'unf', -10, 10);
    nu = grand(1,1,'unf', -10, 10);
    res = besselh1 (nu, arg);
    mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f},\n", "pnl_bessel_h1", "pnl_bessel_h1", arg, nu, real(res), imag(res));
    res = besselh2 (nu, arg);
    mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f, %.18f, %.18f},\n", "pnl_bessel_h2", "pnl_bessel_h2", arg, nu, real(res), imag(res));
end
mclose(dat)

