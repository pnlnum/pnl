
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

dat = mopen ("gamma_test.dat", mode="w");
Lfunc = list( ["gamma", "pnl_sf_gamma"],...
    ["erf", "pnl_sf_erf"], ...
    ["erfc", "pnl_sf_erfc"], ...
    ["gammaln", "pnl_sf_log_gamma"]...
    );

for func=Lfunc
    for i=(1:10)
        arg = grand(1,1,'unf', 0, 10);
        res = evstr(func(1) + "(" + string(arg) + ")"); 
        mfprintf(dat,"{ ""%s"", %s, %.18f, %.18f },\n", func(2), func(2), arg, res)
    end
end

mclose(dat)
