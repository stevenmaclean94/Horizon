fid = fopen('headers.txt');

s = '';
%while (~strcmp(s,'\n'))
    s = [s textscan(fid, '%s')];
%end

tn = {};
tt = {};
nd = [];
c = 1;
for i = 1:size(s{1,1},1)/2
    if s{1,1}{i}(1) ~= 'g'
        tn{c} = s{1,1}{i};
        tt{c} = s{1,1}{i+737};
         nd(:,c) = d(:,i);
       c = c + 1;
    end
end

tn = tn';
tt = tt';
nd = nd';