
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

dat = fopen ("complex_test.dat", mode="w");
Lfunc = ["sqrt", "log", "cos", "sin", "tan", "cosh", "sinh", "tanh"];

for func=Lfunc
    for i=(1:10)
        arg_r = grand(1,1,'uin', 1, 10);
        arg_i = grand(1,1,'uin', 1, 10);
        arg = arg_r + arg_i * %i;
        res = evstr(func + "(" + string(arg_r) + "+ %i *" + string(arg_i) + ")"); 
        fprintf(dat,"{ ""C%s"", C%s, %.18f, %.18f, %.18f, %.18f },\n", func, func, arg_r, arg_i, real(res), imag(res))
    end
end


Lfunc = ["sqrt",  "cos", "sin", "tan", "cosh", "sinh", "tanh"];

for func=Lfunc
    for i=(1:10)
        arg_r = - grand(1,1,'uin', 1, 10);
        arg_i = grand(1,1,'uin', 1, 10);
        arg = arg_r + arg_i * %i;
        res = evstr(func + "(" + string(arg_r) + "+ %i *" + string(arg_i) + ")"); 
        fprintf(dat,"{ ""C%s"", C%s, %.18f, %.18f, %.18f, %.18f },\n", func, func, arg_r, arg_i, real(res), imag(res))
    end
end

Lfunc = ["sqrt",  "cos", "sin", "tan", "cosh", "sinh", "tanh"];

for func=Lfunc
    for i=(1:10)
        arg_r =  grand(1,1,'uin', 1, 10);
        arg_i = - grand(1,1,'uin', 1, 10);
        arg = arg_r + arg_i * %i;
        res = evstr(func + "(" + string(arg_r) + "+ %i *" + string(arg_i) + ")"); 
        fprintf(dat,"{ ""C%s"", C%s, %.18f, %.18f, %.18f, %.18f },\n", func, func, arg_r, arg_i, real(res), imag(res))
    end
end

Lfunc = ["sqrt",  "cos", "sin", "tan", "cosh", "sinh", "tanh"];

for func=Lfunc
    for i=(1:10)
        arg_r = - grand(1,1,'uin', 1, 10);
        arg_i = - grand(1,1,'uin', 1, 10);
        arg = arg_r + arg_i * %i;
        res = evstr(func + "(" + string(arg_r) + "+ %i *" + string(arg_i) + ")"); 
        fprintf(dat,"{ ""C%s"", C%s, %.18f, %.18f, %.18f, %.18f },\n", func, func, arg_r, arg_i, real(res), imag(res))
    end
end
dat.close[];
