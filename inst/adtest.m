function adtest
clear functions;


rand("state",0);
mnpq=[3 4 3 4;4 3 4 3;1 1 1 1;
1 4 3 4;3 1 3 4;3 4 1 4;3 4 3 1;
1 3 4 3;4 1 4 3;4 3 1 3;4 3 4 1;
1 1 3 4;1 1 4 3;3 4 1 1;4 3 1 1];
# plus
for k=1:size(mnpq,1)
  f=@(A,B)A+B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  assert(dCdB,dCdBt,1e-6);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,1e-6);
end
# minus
for k=1:size(mnpq,1)
  f=@(A,B)A-B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  assert(dCdB,dCdBt,1e-6);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,1e-6);
end
## times
for k=1:size(mnpq,1)
  f=@(A,B)A.*B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  assert(dCdB,dCdBt,1e-6);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,1e-6);
end
## rdivide
for k=1:size(mnpq,1)
  f=@(A,B)A./B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
## ldivide
for k=1:size(mnpq,1)
  f=@(A,B)A.\B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
## power
for k=1:size(mnpq,1)
  f=@(A,B)A.^B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
# atan2
for k=1:size(mnpq,1)
  f=@(A,B)atan2(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
# atan2d
for k=1:size(mnpq,1)
  f=@(A,B)atan2d(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
# max
for k=1:size(mnpq,1)
  f=@(A,B)max(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
# min
for k=1:size(mnpq,1)
  f=@(A,B)min(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
# kron
for k=1:size(mnpq,1)
  f=@(A,B)kron(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
## lt
for k=1:size(mnpq,1)
  f=@(A,B)A<B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
## le
for k=1:size(mnpq,1)
  f=@(A,B)A<=B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
## gt
for k=1:size(mnpq,1)
  f=@(A,B)A>B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
## ge
for k=1:size(mnpq,1)
  f=@(A,B)A>=B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
## eq
for k=1:size(mnpq,1)
  f=@(A,B)A==B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
## ne
for k=1:size(mnpq,1)
  f=@(A,B)A~=B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  Ct=f(A,B);
  assert(C,Ct);
  C=f(AA,B);
  assert(C,Ct);
  C=f(A,BB);
  assert(C,Ct);
end
# transpose .'
for k=1:size(mnpq,1)
  f=@(A)A.';
  A=randn(mnpq(k,1),mnpq(k,2));
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
end
# transpose '
for k=1:size(mnpq,1)
  f=@(A)A';
  A=randn(mnpq(k,1),mnpq(k,2));
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
end
# horzcat
for k=1:size(mnpq,1)
  f=@(A)[A A A];
  A=randn(mnpq(k,1),mnpq(k,2));
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
end
# vertcat
for k=1:size(mnpq,1)
  f=@(A)[A;A;A];
  A=randn(mnpq(k,1),mnpq(k,2));
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,1e-6);
end
mnpq=[4 4 4 4;2 4 4 3;
3 4 4 2;1 1 3 4;1 1 4 3;5 3 3 4;
4 3 3 5;4 3 1 1;3 4 1 1;4 4 1 1;
1 1 4 4;1 1 1 1];
for k=1:size(mnpq,1)
  f=@(A,B)A*B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
## mrdivide
mnpq=[4 4 4 4;2 4 3 4;
3 4 2 4;1 1 1 1;4 3 4 3;3 2 4 2];
for k=1:size(mnpq,1)
  f=@(A,B)A/B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
## mldivide
mnpq=[4 4 4 4;4 2 4 3;
4 3 4 2;1 1 1 1;2 4 2 4;3 4 3 1];
for k=1:size(mnpq,1)
  f=@(A,B)A\B;
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-2);
  assert(dCdB,dCdBt,-1e-2);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-2);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-2);
end
## mpower
## Non symetric matrix issues caused by complex normilization of eigvec
## Non positive def different complex values??
## to power of matrix deriv?
mnpq=[2 2 1 1;1 1 1 1];#;1 1 3 3
for k=1:size(mnpq,1)
  f=@(A,B)A^B;
  A=rand(mnpq(k,1),mnpq(k,2));
  B=round(10*rand(mnpq(k,3),mnpq(k,4))-5);
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-1);
  assert(dCdB,dCdBt,-1e-1);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,-1e-1);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-1);
end
## subsref
for s=1
f=@(A,B)A(B);
A=round(10*rand(11,11)+1);
B=round(10*rand(11,11)+1);
AA=ad(A);
Ct=f(A,B);
dCdAt=ad.fd(@(A)f(A,B),A);
C=f(AA,B);
dCdA=C.diff(AA);
assert(C.value,Ct,1e-10);
assert(dCdA,dCdAt,-1e-1);
end
## subsasgn wrong if repeted idx
for s=1
A=1:11;
B=11:-1:1;
idx=struct("type","()");
idx.subs={B};
f=@(A,B)subsasgn(A,idx,B);
AA=ad(A);
BB=ad(B);
C=f(AA,BB);
dCdA=C.diff(AA);
dCdB=C.diff(BB);
Ct=f(A,B);
dCdAt=ad.fd(@(A)f(A,B),A);
dCdBt=ad.fd(@(B)f(A,B),B);
assert(C.value,Ct,1e-10);
assert(dCdA,dCdAt,1e-6);
assert(dCdB,dCdBt,1e-6);
C=f(AA,B);
dCdA=C.diff(AA);
assert(C.value,Ct,1e-10);
assert(dCdA,dCdAt,1e-6);
C=f(A,BB);
dCdB=C.diff(BB);
assert(dCdB,dCdBt,1e-6);
end
## reshape
for s=1
  f=@(A)reshape(A,2,6);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## repmat
for s=1
  f=@(A)repmat(A,2,6);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## repelem
for s=1
  f=@(A)repelem(A,2,6);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## sort
for s=1
  f=@(A)sort(A);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## sortrows
for s=1
  f=@(A)sortrows(A);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## rot90
for s=1
  f=@(A)rot90(A);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## tril
for s=1
  f=@(A)tril(A);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## triu
for s=1
  f=@(A)triu(A);
  A=randn(3,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## vech
for s=1
  f=@(A)vech(A);
  A=randn(4,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## blkdiag
for s=1
  f=@(A)blkdiag(A,A);
  A=randn(4,4);
  Ct=f(A);
  dCdAt=ad.fd(@(A)f(A),A);
  AA=ad(A);
  C=f(AA);
  dCdA=C.diff(AA);
  assert(C.value,Ct,1e-10);
  assert(dCdA,dCdAt,1e-6);
end
## colon
assert(1:3,1:ad(3));
## size
AA=ad(rand(2));
assert(size(AA),size(AA.value));
## end
assert(end(AA,1,1),end(AA.value,1,1));
## columns
assert(columns(AA),columns(AA.value));
## rows
assert(rows(AA),rows(AA.value));
## numel
assert(numel(AA),numel(AA.value));
## ndims
assert(ndims(AA),ndims(AA.value));
## length
assert(length(AA),length(AA.value));
## isempty
assert(isempty(AA),isempty(AA.value));
## isscalar
assert(isscalar(AA),isscalar(AA.value));
## numel
assert(numel(AA),numel(AA.value));
## isnull
assert(isnull(AA),isnull(AA.value));
## isnan
assert(isnan(AA),isnan(AA.value));
## isinf
assert(isinf(AA),isinf(AA.value));
## isnumeric
assert(isnumeric(AA),isnumeric(AA.value));
## single input functions
funcs={"exp"
    "expm1"
    "log"
    "reallog"
    "log10"
    "log2"
    "log1p"
    "sqrt"
    "realsqrt"
    "cbrt"
    "deg2rad"
    "rad2deg"
    "sin"
    "cos"
    "tan"
    "sec"
    "csc"
    "cot"
    "asin"
    "acos"
    "atan"
    "asec"
    "acsc"
    "acot"
    "sinh"
    "cosh"
    "tanh"
    "sech"
    "csch"
    "coth"
    "asinh"
    "acosh"
    "atanh"
    "asech"
    "acsch"
    "acoth"
    "sind"
    "cosd"
    "tand"
    "secd"
    "cscd"
    "cotd"
    "asind"
    "acosd"
    "atand"
    "asecd"
    "acscd"
    "acotd"
    "erf"
    "erfc"
    "erfcx"
    "erfinv"
    "erfcinv"
    "abs"
    "sign"
    "ceil"
    "fix"
    "floor"
    "round"
    "real"
    "imag"
    "conj"
    "arg"
    "max"
    "min"
    "sumsq"
    "mean"
    "median"
    "var"
    "std"
    "sum"
    "cumsum"
    "prod"
    "cumprod"
    "inv"
    "pinv"
    "diag"
    "trace"
    "sqrtm"
    "det"
    "chol"
    "balance"
    "expm"};
AA=ad(randn(2));
for k=1:length(funcs)
  if strcmp(funcs{k},"asec")||strcmp(funcs{k},"asecd")
    AA=ad(rand(2)+1);
  elseif strcmp(funcs{k},"asech")||strcmp(funcs{k},"erfinv")
    AA=ad(rand(2));
  elseif strcmp(funcs{k},"chol")
    A=rand(2);
    A+=A';
    A+=eye(2)*1.1*abs(eigs(A,1,"sa"));
    AA=ad(A);
  elseif strcmp(funcs{k},"reallog")
    A=rand(2);
    AA=ad(A);
  endif
  C=feval(funcs{k},AA);
  Ct=feval(funcs{k},AA.value);
  dCdAt=ad.fd(@(AA)feval(funcs{k},AA),AA.value);
  dCdA=C.diff(AA);
  assert(C.value,Ct,-1e-2);
  assert(dCdA,dCdAt,-1e-2);
end
# dot
mnpq=[2 4 2 4;4 2 4 2];
for k=1:size(mnpq,1)
  f=@(A,B)dot(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
# cross
mnpq=[3 4 3 4;4 3 4 3];
for k=1:size(mnpq,1)
  f=@(A,B)cross(A,B);
  A=randn(mnpq(k,1),mnpq(k,2));
  B=randn(mnpq(k,3),mnpq(k,4));
  AA=ad(A);
  BB=ad(B);
  C=f(AA,BB);
  dCdA=C.diff(AA);
  dCdB=C.diff(BB);
  Ct=f(A,B);
  dCdAt=ad.fd(@(A)f(A,B),A);
  dCdBt=ad.fd(@(B)f(A,B),B);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  assert(dCdB,dCdBt,-1e-4);
  C=f(AA,B);
  dCdA=C.diff(AA);
  assert(C.value,Ct);
  assert(dCdA,dCdAt,-1e-4);
  C=f(A,BB);
  dCdB=C.diff(BB);
  assert(dCdB,dCdBt,-1e-4);
end
#svd
mn=[2 2;3 2;2 3;1 1];
for k=1:size(mn,1)
  AA=ad(randn(mn(k,:)));
  [Ut St Vt]=svd(AA.value,"econ");
  dUdAt=ad.fd(@(AA)nthargout(1:3,@svd,AA,"econ"){1},AA.value);
  dSdAt=ad.fd(@(AA)nthargout(1:3,@svd,AA,"econ"){2},AA.value);
  dVdAt=ad.fd(@(AA)nthargout(1:3,@svd,AA,"econ"){3},AA.value);
  [U S V]=svd(AA);
  assert(U.value,Ut,-1e-2);
  assert(S.value,St,-1e-2);
  assert(V.value,Vt,-1e-2);
  assert(U.diff(AA),dUdAt,-1e-2);
  assert(S.diff(AA),dSdAt,-1e-2);
  assert(V.diff(AA),dVdAt,-1e-2);
end
#qr
for k=1:size(mn,1)
  AA=ad(randn(mn(k,:)));
  [Qt Rt]=qr(AA.value);
  dQdAt=ad.fd(@(AA)nthargout(1:2,@qr,AA){1},AA.value);
  dRdAt=ad.fd(@(AA)nthargout(1:2,@qr,AA){2},AA.value);
  [Q R]=qr(AA);
  assert(Q.value,Qt,-1e-2);
  assert(R.value,Rt,-1e-2);
  assert(Q.diff(AA),dQdAt,-1e-2);
  assert(R.diff(AA),dRdAt,-1e-2);
end
mn=[2 2;1 1];
#eig
for k=1:size(mn,1)
  A=rand(mn(k,:));
  AA=ad(A);
  [Ut St Vt]=eig(AA.value);
  dUdAt=ad.fd(@(AA)nthargout(1:3,@eig,AA){1},AA.value);
  dSdAt=ad.fd(@(AA)nthargout(1:3,@eig,AA){2},AA.value);
  dVdAt=ad.fd(@(AA)nthargout(1:3,@eig,AA){3},AA.value);
  [U S V]=eig(AA);
  assert(U.value,Ut,-1e-2);
  assert(S.value,St,-1e-2);
  assert(V.value,Vt,-1e-2);
  assert(U.diff(AA),dUdAt,-1e-2);
  assert(S.diff(AA),dSdAt,-1e-2);
  assert(V.diff(AA),dVdAt,-1e-2);
end
#eigs
mn=[2 2];
for k=1:size(mn,1)
  A=rand(mn(k,:));
  AA=ad(A);
  [Ut St Vt]=eigs(AA.value,1,"sm");
  dUdAt=ad.fd(@(AA)nthargout(1:3,@eigs,AA,1,"sm"){1},AA.value);
  dSdAt=ad.fd(@(AA)nthargout(1:3,@eigs,AA,1,"sm"){2},AA.value);
  [U S V]=eigs(AA,1,"sm");
  assert(U.value,Ut,-1e-2);
  assert(S.value,St,-1e-2);
  assert(V,Vt,-1e-2);
  assert(U.diff(AA),dUdAt,-1e-2);
  assert(S.diff(AA),dSdAt,-1e-2);
end
#lu
mn=[2 2];
for k=1:size(mn,1)
  A=randn(mn(k,:));
  AA=ad(A);
  [Ut St Vt]=lu(AA.value);
  dUdAt=ad.fd(@(AA)nthargout(1:3,@lu,AA){1},AA.value);
  dSdAt=ad.fd(@(AA)nthargout(1:3,@lu,AA){2},AA.value);
  [U S V]=lu(AA);
  assert(U.value,Ut,-1e-2);
  assert(S.value,St,-1e-2);
  assert(V,Vt,-1e-2);
  assert(U.diff(AA),dUdAt,-1e-2);
  assert(S.diff(AA),dSdAt,-1e-2);
end
## sylvester
A=ad(randn(2));
B=ad(randn(2));
C=ad(randn(2));
D=sylvester(A,B,C);
Dt=sylvester(A.value,B.value,C.value);
dDdAt=ad.fd(@(A)sylvester(A,B.value,C.value),A.value);
dDdBt=ad.fd(@(B)sylvester(A.value,B,C.value),B.value);
dDdCt=ad.fd(@(C)sylvester(A.value,B.value,C),C.value);
assert(D.value,Dt,-1e-2);
assert(D.diff(A),dDdAt,1e-4);
assert(D.diff(B),dDdBt,1e-4);
assert(D.diff(C),dDdCt,1e-4);
## schur
##A=ad(rand(3));
A=ad([0 2 2;0 1 2;1 0 1]);
[Ut St]=schur(A.value);
[U S]=schur(A);
dUdAt=ad.fd(@(A)nthargout(1:2,@schur,A){1},A.value);
dSdAt=ad.fd(@(A)nthargout(1:2,@schur,A){2},A.value);
assert(U.value,Ut,-1e-2);
assert(S.value,St,-1e-2);
assert(U.diff(A),dUdAt,-1e-2);
assert(S.diff(A),dSdAt,-1e-2);
##rsf2csf
[Ur Sr]=rsf2csf(U,S);
[Urt Srt]=rsf2csf(Ut,St);
D=Ur*Sr*Ur';
Dt=Urt*Srt*Urt';
E=Ur'*Ur;
Et=Urt'*Urt;
assert(D.value,Dt,1e-6);
assert(E.value,Et,1e-6);
##norm
A=rand(3);
AA=ad(A);
C=norm(AA);
Ct=norm(A);
dCdAt=ad.fd(@(A)norm(A),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,1);
Ct=norm(A,1);
dCdAt=ad.fd(@(A)norm(A,1),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,inf);
Ct=norm(A,inf);
dCdAt=ad.fd(@(A)norm(A,inf),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,"fro");
Ct=norm(A,"fro");
dCdAt=ad.fd(@(A)norm(A,"fro"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,3);
Ct=norm(A,3);
dCdAt=ad.fd(@(A)norm(A,3),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);

C=norm(AA,2,"rows");
Ct=norm(A,2,"rows");
dCdAt=ad.fd(@(A)norm(A,2,"rows"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,1,"rows");
Ct=norm(A,1,"rows");
dCdAt=ad.fd(@(A)norm(A,1,"rows"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,inf,"rows");
Ct=norm(A,inf,"rows");
dCdAt=ad.fd(@(A)norm(A,inf,"rows"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,"fro","rows");
Ct=norm(A,"fro","rows");
dCdAt=ad.fd(@(A)norm(A,"fro","rows"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,3,"rows");
Ct=norm(A,3,"rows");
dCdAt=ad.fd(@(A)norm(A,3,"rows"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);

p=ad(3);
C=norm(AA,p);
Ct=norm(A,p.value);
dCdAt=ad.fd(@(A)norm(A,p.value),A);
dCdpt=ad.fd(@(p)norm(A,p),p.value);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
assert(C.diff(p),dCdpt,1e-6);
C=norm(AA,p,"rows");
Ct=norm(A,p.value,"rows");
dCdAt=ad.fd(@(A)norm(A,p.value,"rows"),A);
dCdpt=ad.fd(@(p)norm(A,p,"rows"),p.value);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
assert(C.diff(p),dCdpt,1e-6);

C=norm(AA,2,"columns");
Ct=norm(A,2,"columns");
dCdAt=ad.fd(@(A)norm(A,2,"columns"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,1,"columns");
Ct=norm(A,1,"columns");
dCdAt=ad.fd(@(A)norm(A,1,"columns"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,inf,"columns");
Ct=norm(A,inf,"columns");
dCdAt=ad.fd(@(A)norm(A,inf,"columns"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,"fro","columns");
Ct=norm(A,"fro","columns");
dCdAt=ad.fd(@(A)norm(A,"fro","columns"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
C=norm(AA,3,"columns");
Ct=norm(A,3,"columns");
dCdAt=ad.fd(@(A)norm(A,3,"columns"),A);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);

p=ad(3);
C=norm(AA,p,"columns");
Ct=norm(A,p.value,"columns");
dCdAt=ad.fd(@(A)norm(A,p.value,"columns"),A);
dCdpt=ad.fd(@(p)norm(A,p,"columns"),p.value);
assert(C.value,Ct,1e-6);
assert(C.diff(AA),dCdAt,1e-6);
assert(C.diff(p),dCdpt,1e-6);
endfunction
