
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

dat = fopen ("mathtools_test.dat", mode="w");
fprintf(dat,"{ ""pnl_pow_i"", pnl_pow_i, 0., 0, 1 },\n");
fprintf(dat,"{ ""pnl_pow_i"", pnl_pow_i, %f, 0, 1 },\n", exp(1));
fprintf(dat,"{ ""pnl_pow_i"", pnl_pow_i, %f, 0, 1 },\n", -exp(1));

for i=(1:10)
    x = grand(1,1,'nor', 0, 1);
    y = grand(1,1,'uin', 1, 10);
    res = x**y; 
    fprintf(dat,"{ ""pnl_pow_i"", pnl_pow_i, %.18f, %d, %.18f },\n",  x, y, res);
end

for i=(1:10)
    x = grand(1,1,'nor', 0, 1);
    y = -grand(1,1,'uin', 1, 10);
    res = x**y;
    fprintf(dat,"{ ""pnl_pow_i"", pnl_pow_i, %.18f, %d, %.18f },\n",  x, y, res);
end

dat.close[];
