% G1  t1=1; t2=t2+1; t3=t3+1;
disp('G_{1,b+1,c+1}')
clear
syms t t1 t2 t3 b c x;
n=5;

% 商矩阵
t1=1;
t2=b+1;
t3=c+1;

e=2*t1+t2+t3+1;

L=[  2,  0,  0,       1,       1;
     0,  1,  0,       1,       0;
     0,  0,  1,       0,       1;
    t1, t2,  0, t1+t2+1,       1;
    t1,  0, t3,       1, t1+t3+1];

L=eye(n)+t*L;

% 字典序
Z={[1,2],[1,3],[1,4],[1,5],[2,3],[2,4],[2,5],[3,4],[3,5],[4,5]};
%生成2级复合矩阵
m=nchoosek(n,2);
for i=1:m
    for j=1:m
    C(i,j)=det(L(Z{i},Z{j}));
    end
end
%取t的一次项系数，生成2级加性复合矩阵
for i=1:m
    for j=1:m
    dCx = diff(C(i,j),t);
    C(i,j)=subs(dCx,{t},{0});
    end
end

C;
F1=det((x)*eye(m)-C);
F=det((x+e+3-(1.5/(t2+t3+3)))*eye(m)-C);


%F=F*(1024*(b+c+5)^10)
cx=coeffs(F,x);


fileID = fopen('G_{1,b+1,c+1}.txt','w');
fprintf(fileID,'The coefficients of $P(Delta_2(Q^pi(G_{1,b+1,c+1})),x+e(G)+3-1.5/n})$ \n\n');

for i=1:size(cx,2)

    coe=prod(factor(cx(1,i)))*(b + c + 5)^(11-i)*1024;   %保留整系数分子
    fprintf(fileID,['The coefficient of the x^',num2str(i-1),'(numerator only) is:\n' ]);
    fprintf(fileID,'%s\n',coe);

    pos=length(find(coeffs(coe,"all")>0));
    neg=length(find(coeffs(coe,"all")<0));
    fprintf(fileID,['The number of positive and negative terms is respectively:  ',num2str(pos),'   ',num2str(neg),'\n']);
    fprintf(fileID,'\n');
    
end

fclose(fileID);

% c0=cx(1,11)
% c0=factor(c0)
% c0=prod(c0)
% c10=prod(factor(cx(1,11)))
