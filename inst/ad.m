classdef ad<handle
  properties (SetAccess=protected,GetAccess=public)
    value;
    wr2;
    dif;
    id;
    opfunc;
  endproperties

  methods (Static)
    function dcda=diff(c,a)
      if isa(a,"ad")
        if isfield(c.dif,a.id)
          dcda=c.dif.(a.id);
        else
          dcda=[];
        end
      else
        if isfield(c.dif,a)
          dcda=c.dif.(a);
        else
          dcda=[];
        end
      end
    endfunction
    %helpful functions
    function c=boxproduct(a,b)
      [m1 n1]=size(a);
      [m2 n2]=size(b);
##      c=repmat(repelem(a,m2,1),1,n2).*repelem(repmat(b,m1,1),1,n1);#same as below
      c=kron(ones(1,n2),kron(a,ones(m2,1))).*kron(kron(ones(m1,1),b),ones(1,n1));
    end
##    track sparsity
    function [da row col]=difreshape(sizea,numela,numelc,f)
      aidx=reshape(1:numela,sizea(2),sizea(1)).';
      aidx=f(aidx);
      aidx=aidx.'(:);
      aeq=find(aidx.'(:)~=0);
      row=aeq(:);
      aidx(aidx==0)=[];
      col=aidx(:);
      da=sparse(numelc,numela);
      da(sub2ind([numelc,numela],row,col))=1;
    endfunction
##    finite difference
    function dfdx=fd(f,x)
      fx=f(x);
      nvar=numel(x);
      dfdx=zeros(numel(fx),nvar);
      h=eps(max(max(abs(x)(:)),1))^.5;
      for k=1:nvar
        xh=x.';
        xh(k)+=h;
        fxh=(f(xh.')-fx)/h;
        dfdx(:,k)=fxh.'(:);
      end
    endfunction
  endmethods

  methods (Access = protected)
    function c=insertdifn(c,ab)
      wrt={};
      for k=1:length(ab)
        wrt=[wrt ab{k}{1}.wr2];
      end
      c.wr2=unique(wrt);
      for k=1:length(c.wr2)
        di=c.wr2{k};
        c.dif.(di)=sparse(0);
        for j=1:length(ab)
          dadjk=diff(ab{j}{1},di);
          if ~isempty(dadjk)
            c.dif.(di)+=ab{j}{2}*dadjk;
          end
        end
      end
    endfunction
    function c=insertdif2(c,a,b,dcda,dcdb)
      c.wr2=unique([a.wr2 b.wr2]);
      for k=1:length(c.wr2)
        di=c.wr2{k};
        dadak=diff(a,di);
        dbdak=diff(b,di);
        if ~isempty(dadak)&&~isempty(dbdak)
          c.dif.(di)=dcda*dadak+dcdb*dbdak;
        elseif ~isempty(dadak)
          c.dif.(di)=dcda*dadak;
        elseif ~isempty(dbdak)
          c.dif.(di)=dcdb*dbdak;
        end
        dadbk=diff(a,di);
        dbdbk=diff(b,di);
        if ~isempty(dadbk)&&~isempty(dbdbk)
          c.dif.(di)=dcda*dadbk+dcdb*dbdbk;
        elseif ~isempty(dadbk)
          c.dif.(di)=dcda*dadbk;
        elseif ~isempty(dbdbk)
          c.dif.(di)=dcdb*dbdbk;
        end
      end
    endfunction
    function c=insertdif1(c,a,dcda)
      c.wr2=a.wr2;
      for k=1:length(c.wr2)
        di=c.wr2{k};
        dadak=diff(a,di);
        if ~isempty(dadak)
          c.dif.(di)=dcda*dadak;
        end
      end
    endfunction
  endmethods

  methods
    function c=ad(v)
      if isa(v,"ad")
        c.value=v.value;
        c.id=numtimescalled();
        c.dif=v.dif;
        c.dif.(c.id)=speye(numel(v.value));
        c.wr2=unique([v.wr2 {c.id}]);
      else
        c.value=v;
        c.id=numtimescalled();
        c.dif.(c.id)=speye(numel(v));
        c.wr2={c.id};
      end
    endfunction

    %scalar operators
    function c=plus(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(a.value+b.value);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=speye(numela);
          dcdb=speye(numelb);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            if all(~sizeab)
              dcdb=sparse(ones(numela,1));
            elseif sizeab(1)
              dcdb=sparse(repelem(eye(numelb),sizea(2),1));
            else
              dcdb=repmat(speye(numelb),sizea(1),1);
            end
            dcda=speye(numela);
          else
            if all(~sizeab)
              dcda=sparse(ones(numelb,1));
            elseif sizeab(1)
              dcda=sparse(repelem(eye(numela),sizeb(2),1));
            else
              dcda=repmat(speye(numela),sizeb(1),1);
            end
            dcdb=speye(numelb);
          end
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(a.value+b);
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=speye(numela);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            dcda=speye(numela);
          else
            if all(~sizeab)
              dcda=sparse(ones(numelb,1));
            elseif sizeab(1)
              dcda=sparse(repelem(eye(numela),sizeb(2),1));
            else
              dcda=repmat(speye(numela),sizeb(1),1);
            end
          end
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(a+b.value);
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcdb=speye(numelb);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            if all(~sizeab)
              dcdb=sparse(ones(numela,1));
            elseif sizeab(1)
              dcdb=sparse(repelem(eye(numelb),sizea(2),1));
            else
              dcdb=repmat(speye(numelb),sizea(1),1);
            end
          else
            dcdb=speye(numelb);
          end
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="plus";
    endfunction

    function c=uplus(c)
    endfunction

    function c=minus(a,b)
      c=plus(a,uminus(b));
    endfunction

    function c=uminus(a)
      c=ad(-a.value);
      dcda=-speye(numel(a.value));
      c=insertdif1(c,a,dcda);
      c.opfunc="uminus";
    endfunction

    function c=times(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(a.value.*b.value);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=sparse(diag(b.value.'(:)));
          dcdb=sparse(diag(a.value.'(:)));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse(b.value.'(:));
              dcdb=sparse(diag(repmat(a.value,numelb,1)));
            else
              dcda=sparse(diag(repmat(b.value,numela,1)));
              dcdb=sparse(a.value.'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag(repelem(b.value.'(:),sizea(2),1)));
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=a.value.'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=b.value.'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag(repelem(a.value.'(:),sizeb(2),1)));
            end
          else
            if nfout==numela
              dcda=sparse(diag(repmat(b.value.'(:),sizea(1),1)));
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=a.value(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=b.value(:);
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag(repmat(a.value.'(:),sizeb(1),1)));
            end
          end
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(a.value.*b);
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=diag(b.'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse(b.'(:));
            else
              dcda=sparse(diag(repmat(b,numela,1)));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag(repelem(b.'(:),sizea(2),1)));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=b.'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          else
            if nfout==numela
              dcda=sparse(diag(repmat(b.'(:),sizea(1),1)));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=b(:);
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          end
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(a.*b.value);
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcdb=diag(a.'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcdb=sparse(diag(repmat(a,numelb,1)));
            else
              dcdb=sparse(a.'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=a.'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag(repelem(a.'(:),sizeb(2),1)));
            end
          else
            if nfout==numela
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=a(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag(repmat(a.'(:),sizeb(1),1)));
            end
          end
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="times";
    endfunction

    function c=rdivide(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(a.value./b.value);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=diag((1./b.value).'(:));
          dcdb=diag((-a.value./b.value.^2).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse(1./b.value.'(:));
              dcdb=sparse(diag((-a.value./b.value.^2).'(:)));
            else
              dcda=sparse(diag(repmat(1./b.value,numela,1)));
              dcdb=sparse((-a.value./b.value.^2).'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag(repelem(1./b.value.'(:),sizea(2),1)));
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(-a.value./b.value.^2).'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=1./b.value.'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag((-a.value./b.value.^2).'(:)));
            end
          else
            if nfout==numela
              dcda=sparse(diag(repmat(1./b.value.'(:),sizea(1),1)));
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=(-a.value./b.value.^2)(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=1./b.value(:);
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag((-a.value./b.value.^2).'(:)));
            end
          end
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(a.value./b);
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=diag((1./b).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse(1./b.'(:));
            else
              dcda=sparse(diag(repmat(1./b,numela,1)));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag(repelem(1./b.'(:),sizea(2),1)));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=1./b.'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          else
            if nfout==numela
              dcda=sparse(diag(repmat(1./b.'(:),sizea(1),1)));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=1./b(:);
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          end
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(a./b.value);
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcdb=diag((-a./b.value.^2).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcdb=sparse(diag((-a./b.value.^2).'(:)));
            else
              dcdb=sparse((-a./b.value.^2).'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(-a./b.value.^2).'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag((-a./b.value.^2).'(:)));
            end
          else
            if nfout==numela
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=(-a./b.value.^2)(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag((-a./b.value.^2).'(:)));
            end
          end
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="rdivide";
    endfunction

    function c=ldivide(a,b)
      c=rdivide(b,a);
      c.opfunc="ldivide";
    endfunction

    function c=power(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(a.value.^b.value);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=diag((b.value.*a.value.^(b.value-1)).'(:));
          dcdb=diag((log(a.value).*a.value.^b.value).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse((b.value.*a.value.^(b.value-1)).'(:));
              dcdb=sparse(diag((log(a.value).*a.value.^b.value).'(:)));
            else
              dcda=sparse(diag((b.value.*a.value.^(b.value-1)).'(:)));
              dcdb=sparse((log(a.value).*a.value.^b.value).'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag((b.value.*a.value.^(b.value-1)).'(:)));
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(log(a.value).*a.value.^b.value).'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=(b.value.*a.value.^(b.value-1)).'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag((log(a.value).*a.value.^b.value).'(:)));
            end
          else
            if nfout==numela
              dcda=sparse(diag((b.value.*a.value.^(b.value-1)).'(:)));
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=(log(a.value).*a.value.^b.value)(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=(b.value.*a.value.^(b.value-1));
              dcda=sparse(reshape(dcda,numelb,numela));
              dcdb=sparse(diag((log(a.value).*a.value.^b.value).'(:)));
            end
          end
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(a.value.^b);
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=diag((b.*a.value.^(b-1)).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcda=sparse((b.*a.value.^(b-1)).'(:));
            else
              dcda=sparse(diag((b.*a.value.^(b-1)).'(:)));
            end
          elseif sizeab(1)
            if nfout==numela
              dcda=sparse(diag((b.*a.value.^(b-1)).'(:)));
            else
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=(b.*a.value.^(b-1)).'(:);
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          else
            if nfout==numela
              dcda=sparse(diag((b.*a.value.^(b-1)).'(:)));
            else
              dcda(logical(repmat(eye(numela),sizeb(1),1)))=(b.*a.value.^(b-1));
              dcda=sparse(reshape(dcda,numelb,numela));
            end
          end
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(a.^b.value);
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcdb=diag((log(a).*a.^b.value).'(:));
        else
          nfout=max(numela,numelb);
          if all(~sizeab)
            if numela==1
              dcdb=sparse(diag((log(a).*a.^b.value).'(:)));
            else
              dcdb=sparse((log(a).*a.^b.value).'(:));
            end
          elseif sizeab(1)
            if nfout==numela
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(log(a).*a.^b.value).'(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag((log(a).*a.^b.value).'(:)));
            end
          else
            if nfout==numela
              dcdb(logical(repmat(eye(numelb),sizea(1),1)))=(log(a).*a.^b.value)(:);
              dcdb=sparse(reshape(dcdb,numela,numelb));
            else
              dcdb=sparse(diag((log(a).*a.^b.value).'(:)));
            end
          end
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="power";
    endfunction

    %Matrix operators
    function c=mtimes(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(a.value*b.value);
        sizea=size(a.value);
        sizeb=size(b.value);
        if all(sizea==1)||all(sizeb==1)
          if all(sizea==1)
            dcda=sparse(b.value.'(:));
            dcdb=a.value*speye(prod(sizeb));
          else
            dcda=b.value*speye(prod(sizea));
            dcdb=sparse(a.value.'(:));
          end
        else
          dcda=kron(speye(size(a.value,1)),b.value.');
          dcdb=kron(a.value,speye(size(b.value,2)));
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(a.value*b);
        sizea=size(a.value);
        sizeb=size(b);
        if all(sizea==1)||all(sizeb==1)
          if all(sizea==1)
            dcda=sparse(b.'(:));
          else
            dcda=b*speye(prod(sizea));
          end
        else
          dcda=kron(speye(size(a.value,1)),b.');
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(a*b.value);
        sizea=size(a);
        sizeb=size(b.value);
        if all(sizea==1)||all(sizeb==1)
          if all(sizea==1)
            dcdb=a*speye(prod(sizeb));
          else
            dcdb=sparse(a.'(:));
          end
        else
          dcdb=kron(a,speye(size(b.value,2)));
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="mtimes";
    endfunction

    function c=mrdivide(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        if all(size(b.value)==1)
          c=ad(a.value/b.value);
          dcda=speye(numel(a))/b.value;
          dcdb=sparse(-a.value.'(:)/b.value^2);
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif size(b.value,1)==size(b.value,2)
          lastwarn("");
          c=ad(a.value/b.value);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            invb=pinv(b.value);
          else
            invb=inv(b.value);
          end
          dcda=kron(speye(size(a.value,1)),invb.');
          dcdb=kron(a.value*-invb,sparse(invb).');
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif size(b.value,1)<size(b.value,2)
          c=(a*b')/(b*b');
        else
          c=a*pinv(b);
        end
      elseif isa(a,"ad")
        if all(size(b)==1)
          c=ad(a.value/b);
          dcda=speye(numel(a))/b;
          c=insertdif1(c,a,dcda);
        elseif size(b,1)==size(b,2)
          lastwarn("");
          c=ad(a.value/b);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            invb=pinv(b);
          else
            invb=inv(b);
          end
          dcda=kron(speye(size(a.value,1)),invb.');
          c=insertdif1(c,a,dcda);
        elseif size(b,1)<size(b,2)
          c=(a*b')/(b*b');
        else
          c=a*pinv(b);
        end
      elseif isa(b,"ad")
        if all(size(b.value)==1)
          c=ad(a/b.value);
          dcdb=sparse(-a.'(:)/b.value^2);
          c=insertdif1(c,b,dcdb);
        elseif size(b.value,1)==size(b.value,2)
          lastwarn("");
          c=ad(a/b.value);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            invb=pinv(b.value);
          else
            invb=inv(b.value);
          end
          dcdb=kron(a*-invb,sparse(invb).');
          c=insertdif1(c,b,dcdb);
        elseif size(b.value,1)<size(b.value,2)
          c=(a*b')/(b*b');
        else
          c=a*pinv(b);
        end
      end
      c.opfunc="mrdivide";
    endfunction

    function c=mldivide(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        if all(size(a.value)==1)
          c=ad(a.value\b.value);
          dcda=sparse(-a.value^2\b.value.'(:));
          dcdb=a.value\speye(numel(b.value));
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif size(a.value,1)==size(a.value,2)
          lastwarn("");
          c=ad(a.value\b.value);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            inva=pinv(a.value);
          else
            inva=inv(a.value);
          end
          dcda=kron(-inva,sparse(inva*b.value).');
          dcdb=kron(inva,speye(size(b.value,2)));
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif size(a.value,1)>size(a.value,2)
          c=((a'*a)\(a'*b));
        else
          c=pinv(a)*b;
        end
      elseif isa(a,"ad")
        if all(size(a.value)==1)
          c=ad(a.value\b);
          dcda=sparse(-a.value^2\b.'(:));
          c=insertdif1(c,a,dcda);
        elseif size(a.value,1)==size(a.value,2)
          lastwarn("");
          c=ad(a.value\b);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            inva=pinv(a.value);
          else
            inva=inv(a.value);
          end
          dcda=kron(-inva,sparse(inva*b).');
          c=insertdif1(c,a,dcda);
        elseif size(a.value,1)>size(a.value,2)
          c=((a'*a)\(a'*b));
        else
          c=pinv(a)*b;
        end
      elseif isa(b,"ad")
        if all(size(a)==1)
          c=ad(a\b.value);
          dcdb=a\speye(numel(b.value));
          c=insertdif1(c,b,dcdb);
        elseif size(a,1)==size(a,2)
          lastwarn("");
          c=ad(a\b.value);
          [~,msgid]=lastwarn();
          if strcmp(msgid,"Octave:singular-matrix")||strcmp(msgid,"Octave:nearly-singular-matrix")
            inva=pinv(a);
          else
            inva=inv(a);
          end
          dcdb=kron(inva,speye(size(b.value,2)));
          c=insertdif1(c,b,dcdb);
        elseif size(a,1)>size(a,2)
          c=((a'*a)\(a'*b));
        else
          c=pinv(a)*b;
        end
      end
      c.opfunc="mldivide";
    endfunction

    function c=mpower(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        if all(size(a.value)==1)&&all(size(b.value)==1)
          c=ad(a.value^b.value);
          dcda=b.value*a.value^(b.value-1);
          dcdb=log(a.value)*a.value^b.value;
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif all(size(b.value)==1)&&rem(b.value*2,1)==0
          c=ad(a.value^b.value);
          if rem(b.value,1)==0
            bb=abs(b.value);
          else
            bb=abs(b.value)*2;
          end
          if b.value<0
            aa=inv(a.value);
          else
            aa=a.value;
          end
          difa=zeros(numel(aa));
          for i=0:bb-1
            difa+=kron(aa^i,(aa^(bb-1-i)).');
          end
          if b.value<0
            difa=difa*kron(-aa,aa.');
          end
          if ~rem(b.value,1)==0
            aab=aa^abs(b.value);
            difa=difa/(kron(eye(size(aa)),aab.')+kron(aab,eye(size(aa))));
          end
          dcda=sparse(difa);
          dcdb=sparse(logm(a.value)*a.value^b.value).'(:);
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif all(size(b.value)==1)
          [V L]=eig(a);
          diagL=diag(L);
          diagL=diagL.^b;
          L=diag(diagL);
          c=V*L/V;
        else
          [V L]=eig(b);
          diagL=diag(L);
          diagL=a.^diagL;
          L=diag(diagL);
          c=V*L/V;
        end
      elseif isa(a,"ad")
        c=ad(a.value^b);
        if all(size(a.value)==1)&&all(size(b)==1)
          dcda=b*a.value^(b-1);
          c=insertdif1(c,a,dcda);
        elseif all(size(b)==1)&&rem(b*2,1)==0
          if rem(b,1)==0
            bb=abs(b);
          else
            bb=abs(b)*2;
          end
          if b<0
            aa=inv(a.value);
          else
            aa=a.value;
          end
          difa=zeros(numel(aa));
          for i=0:bb-1
            difa+=kron(aa^i,(aa^(bb-1-i)).');
          end
          if b<0
            difa=difa*kron(-aa,aa.');
          end
          if ~rem(b,1)==0
            aab=aa^abs(b);
            difa=difa/(kron(eye(size(aa)),aab.')+kron(aab,eye(size(aa))));
          end
          dcda=sparse(difa);
          c=insertdif1(c,a,dcda);
        elseif all(size(b)==1)
          [V L]=eig(a);
          diagL=diag(L);
          diagL=diagL.^b;
          L=diag(diagL);
          c=V*L/V;
        else
          [V L]=eig(b);
          diagL=diag(L);
          diagL=a.^diagL;
          L=diag(diagL);
          c=V*L/V;
        end
      elseif isa(b,"ad")
        c=ad(a^b.value);
        if all(size(b.value)==1)
          dcdb=sparse(logm(a)*a^b.value).'(:);
          c=insertdif1(c,b,dcdb);
        else
          [V L]=eig(b);
          diagL=diag(L);
          diagL=a.^diagL;
          L=diag(diagL);
          c=V*L/V;
        end
      end
      c.opfunc="mpower";
    endfunction

    %conditions
    function bool=lt(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value<b.value;
      elseif isa(a,"ad")
        bool=a.value<b;
      elseif isa(b,"ad")
        bool=a<b.value;
      end
    endfunction

    function bool=le(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value<=b.value;
      elseif isa(a,"ad")
        bool=a.value<=b;
      elseif isa(b,"ad")
        bool=a<=b.value;
      end
    endfunction

    function bool=gt(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value>b.value;
      elseif isa(a,"ad")
        bool=a.value>b;
      elseif isa(b,"ad")
        bool=a>b.value;
      end
    endfunction

    function bool=ge(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value>=b.value;
      elseif isa(a,"ad")
        bool=a.value>=b;
      elseif isa(b,"ad")
        bool=a>=b.value;
      end
    endfunction

    function bool=eq(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value==b.value;
      elseif isa(a,"ad")
        bool=a.value==b;
      elseif isa(b,"ad")
        bool=a==b.value;
      end
    endfunction

    function bool=ne(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        bool=a.value~=b.value;
      elseif isa(a,"ad")
        bool=a.value~=b;
      elseif isa(b,"ad")
        bool=a~=b.value;
      end
    endfunction

    function bool=isdiag(a)
      bool=isdiag(a.value);
    endfunction

    %reshaping, reordering, concentrating and subsituting values
    function c=ctranspose(a)
      c=ad(a.value');
      dcda=sparse(ad.boxproduct(eye(size(a.value,2)),eye(size(a.value,1))));
      c=insertdif1(c,a,dcda);
      c.opfunc="ctranspose";
    endfunction

    function c=transpose(a)
      c=ad(a.value.');
      dcda=sparse(ad.boxproduct(eye(size(a.value,2)),eye(size(a.value,1))));
      c=insertdif1(c,a,dcda);
      c.opfunc="transpose";
    endfunction

    function c=horzcat(a,b,varargin)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad([a.value b.value]);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        numelc=numel(c.value);
        dcda=ad.difreshape(sizea,numela,numelc,@(x)[x sparse(sizeb(1),sizeb(2))]);
        dcdb=ad.difreshape(sizeb,numelb,numelc,@(x)[sparse(sizea(1),sizea(2)) x]);
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad([a.value b]);
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        numelc=numel(c.value);
        dcda=ad.difreshape(sizea,numela,numelc,@(x)[x sparse(sizeb(1),sizeb(2))]);
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad([a b.value]);
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        numelc=numel(c.value);
        dcdb=ad.difreshape(sizeb,numelb,numelc,@(x)[sparse(sizea(1),sizea(2)) x]);
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="horzcat";
      if length(varargin)>0
        c=horzcat(c,varargin{1},varargin{2:end});
        c.opfunc="horzcat";
      end
    endfunction

    function c=vertcat(a,b,varargin)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad([a.value;b.value]);
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        dcda=[speye(numela);zeros(numelb,numela)];
        dcdb=[zeros(numela,numelb);speye(numelb)];
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad([a.value;b]);
        numela=numel(a.value);
        numelb=numel(b);
        dcda=[speye(numela);zeros(numelb,numela)];
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad([a;b.value]);
        numela=numel(a);
        numelb=numel(b.value);
        dcdb=[zeros(numela,numelb);speye(numelb)];
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="vertcat";
      if length(varargin)>0
        c=vertcat(c,varargin{1},varargin{2:end});
      end
    endfunction

    function c=subsref(a,s)
      sorig=s;
      s=sorig(1);
      try
        s.type=="()";
      catch
        c=builtin('subsref',a,sorig);
        return
      end
      if s.type=="()"
        c=ad(a.value(s.subs{:}));
        sizea=size(a.value);
        numela=numel(a.value);
        numelc=numel(c.value);
        dcda=ad.difreshape(sizea,numela,numelc,@(x)subsref(x,s));
        c=insertdif1(c,a,dcda);
        c.opfunc="subsref";
        if length(sorig)>1
          c=subsref(c,sorig(2:end));
        end
      elseif s.type=="."
        if strcmp(s.subs,"diff")
          c=ad.diff(a,sorig(2).subs{:});
        else
          c=builtin('subsref',a,sorig);
        end
      else
        c=builtin('subsref',a,sorig);
      end
    endfunction

    function c=subsasgn(a,s,b)
      if s.type=="()"
        if isa(a,"ad")&&isa(b,"ad")
          val=a.value;
          val(s.subs{:})=b.value;
          c=ad(val);
          sizea=size(a.value);
          sizeb=size(b.value);
          sizec=size(c.value);
          numela=numel(a.value);
          numelb=numel(b.value);
          numelc=numel(c.value);
          aidx=reshape(1:numela,sizea(2),sizea(1)).';
          aidx(s.subs{:})=sparse(sizeb(1),sizeb(2));
          aidx=aidx.'(:);
          aeq=find(aidx.'(:)~=0);
          aidx(aidx==0)=[];
          dcda=sparse(numelc,numela);
          dcda(sub2ind([numelc,numela],aeq(:),aidx(:)))=1;
          cidx=reshape(1:numelc,sizec(2),sizec(1)).';
          cidx=cidx(s.subs{:})(:);
          dcdb=sparse(numelc,numelb);
          dcdb(sub2ind([numelc,numelb],cidx,reshape(1:numelb,sizeb(2),sizeb(1)).'(:)))=1;
          c=insertdif2(c,a,b,dcda,dcdb);
        elseif isa(a,"ad")
          val=a.value;
          val(s.subs{:})=b;
          c=ad(val);
          sizea=size(a.value);
          sizeb=size(b);
          sizec=size(c.value);
          numela=numel(a.value);
          numelc=numel(c.value);
          cidx=reshape(1:numelc,sizec(2),sizec(1)).';
          cidx=cidx(s.subs{:})(:);
          aidx=reshape(1:numela,sizea(2),sizea(1)).';
          aidx(s.subs{:})=sparse(sizeb(1),sizeb(2));
          aidx=aidx.'(:);
          aeq=find(aidx.'(:)~=0);
          aidx(aidx==0)=[];
          dcda=sparse(numelc,numela);
          dcda(sub2ind([numelc,numela],aeq(:),aidx(:)))=1;
          c=insertdif1(c,a,dcda);
        elseif isa(b,"ad")
          %does not get called tries to convert object to array: workaround create automatic differentiation object first and join variable or create variable through concentration
          val=a;
          val(s.subs{:})=b.value;
          c=ad(val);
          sizea=size(a);
          sizeb=size(b.value);
          sizec=size(c.value);
          numelb=numel(b.value);
          numelc=numel(c.value);
          cidx=reshape(1:numelc,sizec(2),sizec(1)).';
          cidx=cidx(s.subs{:})(:);
          dcdb=sparse(numelc,numelb);
          dcdb(sub2ind([numelc,numelb],cidx,reshape(1:numelb,sizeb(2),sizeb(1)).'(:)))=1;
          c=insertdif1(c,b,dcdb);
        end
        c.opfunc="subsasgn";
      end
    endfunction

    function c=reshape(a,varargin)
      c=ad(reshape(a.value,varargin{:}));
      sizea=size(a.value);
      sizec=size(c.value);
      numela=numel(a.value);
      numelc=numel(c.value);
      dcda=ad.difreshape(sizea,numela,numelc,@(x)reshape(x,varargin{:}));
      c=insertdif1(c,a,dcda);
      c.opfunc="reshape";
    endfunction

    function c=repmat(a,varargin)
      c=ad(repmat(a.value,varargin{:}));
      sizea=size(a.value);
      sizec=size(c.value);
      numela=numel(a.value);
      numelc=numel(c.value);
      dcda=ad.difreshape(sizea,numela,numelc,@(x)repmat(x,varargin{:}));
      c=insertdif1(c,a,dcda);
      c.opfunc="repmat";
    endfunction

    function c=repelem(a,varargin)
      c=ad(repelem(a.value,varargin{:}));
      sizea=size(a.value);
      sizec=size(c.value);
      numela=numel(a.value);
      numelc=numel(c.value);
      dcda=ad.difreshape(sizea,numela,numelc,@(x)repelem(x,varargin{:}));
      c=insertdif1(c,a,dcda);
      c.opfunc="repelem";
    endfunction

    function [c,I]=sort(a,DIM,MODE)
      sizea=size(a.value);
      if nargin()<2
        MODE="ascend";
      elseif nargin()<3&&isa(DIM,"char")
        MODE=DIM;
      elseif nargin()<3
        MODE="ascend";
      end
      if nargin()<2||~isa(DIM,"numeric")||length(DIM)~=1||isempty(DIM)
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      [~,I]=sort(a.value,DIM,MODE);
      idx.type="()";
      if DIM==1
        idx.subs={I+(0:sizea(1):sizea(1)*(sizea(2)-1))};
      elseif DIM==2
        idx.subs={reshape(sub2ind(size(a.value),repmat((1:sizea(1)).',sizea(2),1),I(:)),size(I))};
      end
      c=subsref(a,idx);
      c.opfunc="sort";
    endfunction

    function [c,I]=sortrows(a,C)
      if nargin()<2
        [~,I]=sortrows(a.value);
      else
        [~,I]=sortrows(a.value,C);
      end
      idx.type="()";
      idx.subs={I,":"};
      c=subsref(a,idx);
      c.opfunc="sortrows";
    endfunction

    function c=rot90(a,varargin)
      c=ad(rot90(a.value,varargin{:}));
      sizea=size(a.value);
      numela=numel(a.value);
      numelc=numel(c.value);
      dcda=ad.difreshape(sizea,numela,numelc,@(x)rot90(x,varargin{:}));
      c=insertdif1(c,a,dcda);;
      c.opfunc="rot90";
    endfunction

    function c=tril(a,K,PACK)
      sizea=size(a.value);
      numela=numel(a.value);
      if nargin()<2
        K=0;
      end
      c=ad(tril(a.value,K));
      d=zeros(numela,1);
      idx=tril(reshape(1:numela,fliplr(sizea)).',K);
      idx(idx==0)=[];
      d(idx)=1;
      dcda=spdiags(d,0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="tril";
      if nargin()==3&&strcmp(PACK,"pack")
        i.type="()";
        i.subs={":"};
        c=subsref(c,i);
      end
    endfunction

    function c=triu(a,K,PACK)
      sizea=size(a.value);
      numela=numel(a.value);
      if nargin()<2
        K=0;
      end
      c=ad(triu(a.value,K));
      d=zeros(numela,1);
      idx=triu(reshape(1:numela,fliplr(sizea)).',K);
      idx(idx==0)=[];
      d(idx)=1;
      dcda=spdiags(d,0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="triu";
      if nargin()==3&&strcmp(PACK,"pack")
        i.type="()";
        i.subs={":"};
        c=subsref(c,i);
      end
    endfunction

    function c=vech(a)
      sizea=size(a.value);
      numela=numel(a.value);
      idx=tril(reshape(1:numela,fliplr(sizea)));
      idx(idx==0)=[];
      i.type="()";
      i.subs={idx};
      c=subsref(a,i);
      c.opfunc="vech";
    endfunction

    function c=blkdiag(varargin)
      ninputs=length(varargin);
      sizes=zeros(ninputs,2);
      numels=zeros(ninputs,1);
      for k=1:ninputs
        sizes(k,:)=size(varargin{k});
        numels(k)=numel(varargin{k});
      end
      csizes=cumsum(sizes,1);
      m=zeros(csizes(end,:));
      numelm=numel(m);
      d={};
      outputline=reshape(1:numelm,fliplr(csizes(end,:))).';
      for k=1:ninputs
        if k==1
          if isa(varargin{k},"ad")
            m(1:csizes(k,1),1:csizes(k,2))=varargin{k}.value;
            d{end+1}={varargin{k} sparse(numelm,numels(k))};
            inputline=reshape(1:numels(k),fliplr(sizes(k,:))).'(:);
            d{end}{2}(sub2ind([numelm,numels(k)],outputline(1:csizes(k,1),1:csizes(k,2))(:),inputline))=1;
          else
            m(1:csizes(k,1),1:csizes(k,2))=varargin{k};
          end
        else
          if isa(varargin{k},"ad")
            m(csizes(k-1,1)+1:csizes(k,1),csizes(k-1,2)+1:csizes(k,2))=varargin{k}.value;
            d{end+1}={varargin{k} sparse(numelm,numels(k))};
            inputline=reshape(1:numels(k),fliplr(sizes(k,:))).'(:);
            d{end}{2}(sub2ind([numelm,numels(k)],outputline(csizes(k-1,1)+1:csizes(k,1),csizes(k-1,2)+1:csizes(k,2))(:),inputline))=1;
          else
            m(csizes(k-1,1)+1:csizes(k,1),csizes(k-1,2)+1:csizes(k,2))=varargin{k};
          end
        end
      end
      c=ad(m);
      c=insertdifn(c,d);
      c.opfunc="blkdiag";
    endfunction

    function r=colon(base,increment,limit)
      if nargin==2
        if isa(base,'ad')
          base=base.value;
        endif
        if isa(increment,'ad')
          limit=increment.value;
        endif
        r=base:limit;
      else
        if isa(base,'ad')
          base=base.value;
        endif
        if isa(increment,'ad')
          increment=increment.value;
        endif
        if isa(limit,'ad')
          limit=limit.value;
        endif
        r=base:increment:limit;
      endif
    endfunction

    %matrix properties
    function varargout=size(a,varargin)
      nout=nargout;
      if nout<=1
        varargout={size(a.value,varargin{:})};
      else
        varargout=nthargout(1:max(nargout,1),@size,a.value,varargin{:});
      end
    end
    function lastindex=end(a,k,n)
      if n==1
        lastindex=numel(a.value);
      else
        lastindex=size(a.value,k);
      endif
    endfunction
    function v=columns(a)
      v=columns(a.value);
    endfunction
    function v=rows(a)
      v=rows(a.value);
    endfunction
    function v=numel(a,varargin)
      v=numel(a.value,varargin{:});
    endfunction
    function v=ndims(a)
      v=ndims(a.value);
    end
    function v=length(a)
      v=length(a.value);
    end
    function v=isempty(a)
      v=isempty(a.value);
    end
    function v=isscalar(a)
      v=all(size(a.value)==1);
    end
    function v=isnull(a)
      v=isnull(a.value);
    end
    function v=isnan(a)
      v=isnan(a.value);
    end
    function v=isinf(a)
      v=isinf(a.value);
    end
    function bool=isnumeric(a)
      bool=isnumeric(a.value);
    end
##    function bool=isreal(a)
##      bool=isreal(a.value);
##    endfunction
##    function v=isvector(a)
##      v=ndims(a.value)==2&&any(size(a.value)==1);
##    end

    %How to display the class
    function display(a)
      printf("%s =\n",inputname(1));
      disp(a.value);
      disp("Automatic Differentiation variable.");
      printf("To obtain the function value type: %s.value\n",inputname(1));
      printf("To obtain the derivative type dfdx=ad.diff(f,x) or dfdx=%s.diff(x)\n",inputname(1));
    endfunction

    %Scalar functions
    function c=exp(a)
      c=ad(exp(a.value));
      numela=numel(a.value);
      dcda=spdiags(c.value.'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="exp";
    endfunction

    function c=expm1(a)
      c=ad(expm1(a.value));
      numela=numel(a.value);
      dcda=spdiags(1+c.value.'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="expm1";
    endfunction

    function c=log(a)
      c=ad(log(a.value));
      numela=numel(a.value);
      dcda=spdiags(1./a.value.'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="log";
    endfunction

    function c=reallog(a)
      c=ad(reallog(a.value));
      numela=numel(a.value);
      dcda=spdiags(1./a.value.'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="reallog";
    endfunction

    function c=log10(a)
      c=ad(log10(a.value));
      numela=numel(a.value);
      dcda=spdiags((1./(log(10)*a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="log10";
    endfunction

    function [c E]=log2(a)
      numela=numel(a.value);
      if nargout()<=1
        c=ad(log2(a.value));
        dcda=spdiags((1./(log(2)*a.value)).'(:),0,numela,numela);
        c=insertdif1(c,a,dcda);
        c.opfunc="log2";
      elseif nargout()==2
        [F E]=log2(a.value);
        c=ad(F);
        dcda=spdiags(1./(2.^E).'(:),0,numela,numela);
        c=insertdif1(c,a,dcda);
        c.opfunc="log2";
      end
    endfunction

    function c=log1p(a)
      numela=numel(a.value);
      c=ad(log1p(a.value));
      dcda=spdiags((1./(1+a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="log1p";
    endfunction

    function c=sqrt(a)
      numela=numel(a.value);
      c=ad(sqrt(a.value));
      dcda=spdiags((.5./sqrt(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sqrt";
    endfunction

    function c=realsqrt(a)
      numela=numel(a.value);
      c=ad(realsqrt(a.value));
      dcda=spdiags((.5./sqrt(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="realsqrt";
    endfunction

    function c=cbrt(a)
      numela=numel(a.value);
      c=ad(cbrt(a.value));
      dcda=spdiags((1./(3*(a.value).^(2/3))).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="realsqrt";
    endfunction

    function c=deg2rad(a)
      numela=numel(a.value);
      c=ad(deg2rad(a.value));
      dcda=spdiags(pi/180*ones(size(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="deg2rad";
    endfunction

    function c=rad2deg(a)
      numela=numel(a.value);
      c=ad(rad2deg(a.value));
      dcda=spdiags(180/pi*ones(size(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="rad2deg";
    endfunction

    function c=sin(a)
      numela=numel(a.value);
      c=ad(sin(a.value));
      dcda=spdiags(cos(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sin";
    endfunction

    function c=cos(a)
      numela=numel(a.value);
      c=ad(cos(a.value));
      dcda=spdiags(-sin(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cos";
    endfunction

    function c=tan(a)
      numela=numel(a.value);
      c=ad(tan(a.value));
      dcda=spdiags((1+tan(a.value).^2).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="tan";
    endfunction

    function c=sec(a)
      numela=numel(a.value);
      c=ad(sec(a.value));
      dcda=spdiags((tan(a.value).*sec(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sec";
    endfunction

    function c=csc(a)
      numela=numel(a.value);
      c=ad(csc(a.value));
      dcda=spdiags((-cot(a.value).*csc(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="csc";
    endfunction

    function c=cot(a)
      numela=numel(a.value);
      c=ad(cot(a.value));
      dcda=spdiags((-cot(a.value).^2-1).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cot";
    endfunction

    function c=asin(a)
      numela=numel(a.value);
      c=ad(asin(a.value));
      dcda=spdiags((1./(1-a.value.^2).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asin";
    endfunction

    function c=acos(a)
      numela=numel(a.value);
      c=ad(acos(a.value));
      dcda=spdiags((-1./(1-a.value.^2).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acos";
    endfunction

    function c=atan(a)
      numela=numel(a.value);
      c=ad(atan(a.value));
      dcda=spdiags((1./(1+a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="atan";
    endfunction

    function c=asec(a)
      numela=numel(a.value);
      c=ad(asec(a.value));
      dcda=spdiags((1./(a.value.^2.*(1-1./a.value.^2).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asec";
    endfunction

    function c=acsc(a)
      numela=numel(a.value);
      c=ad(acsc(a.value));
      dcda=spdiags((-1./(a.value.^2.*(1-1./a.value.^2).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acsc";
    endfunction

    function c=acot(a)
      numela=numel(a.value);
      c=ad(acot(a.value));
      dcda=spdiags((-1./(a.value.^2+1)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acot";
    endfunction

    function c=sinh(a)
      numela=numel(a.value);
      c=ad(sinh(a.value));
      dcda=spdiags(cosh(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sinh";
    endfunction

    function c=cosh(a)
      numela=numel(a.value);
      c=ad(cosh(a.value));
      dcda=spdiags(sinh(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cosh";
    endfunction

    function c=tanh(a)
      numela=numel(a.value);
      c=ad(tanh(a.value));
      dcda=spdiags((1-tanh(a.value).^2).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="tanh";
    endfunction

    function c=sech(a)
      numela=numel(a.value);
      c=ad(sech(a.value));
      dcda=spdiags((-tanh(a.value).*sech(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sech";
    endfunction

    function c=csch(a)
      numela=numel(a.value);
      c=ad(csch(a.value));
      dcda=spdiags((-coth(a.value).*csch(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sech";
    endfunction

    function c=coth(a)
      numela=numel(a.value);
      c=ad(coth(a.value));
      dcda=spdiags((-1./sinh(a.value).^2).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="coth";
    endfunction

    function c=asinh(a)
      numela=numel(a.value);
      c=ad(asinh(a.value));
      dcda=spdiags((1./(1+a.value.^2).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asinh";
    endfunction

    function c=acosh(a)
      numela=numel(a.value);
      c=ad(acosh(a.value));
      dcda=spdiags((1./(a.value.^2-1).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acosh";
    endfunction

    function c=atanh(a)
      numela=numel(a.value);
      c=ad(atanh(a.value));
      dcda=spdiags((1./(1-a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="atanh";
    endfunction

    function c=asech(a)
      numela=numel(a.value);
      c=ad(asech(a.value));
      dcda=spdiags((-1./(a.value.^2.*(1./a.value-1).^.5.*(1./a.value+1).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asech";
    endfunction

    function c=acsch(a)
      numela=numel(a.value);
      c=ad(acsch(a.value));
      dcda=spdiags((-1./(a.value.^2.*(1+1./a.value.^2).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acsch";
    endfunction

    function c=acoth(a)
      numela=numel(a.value);
      c=ad(acoth(a.value));
      dcda=spdiags((1./(1-a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acoth";
    endfunction

    function c=atan2(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(atan2(a.value,b.value));
        sizea=size(a.value);
        sizeb=size(b.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        numelc=numel(c.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=spdiags((b.value./(a.value.^2+b.value.^2)).'(:),0,numelc,numelc);
          dcdb=spdiags((-a.value./(a.value.^2+b.value.^2)).'(:),0,numelc,numelc);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            if all(~sizeab)
              dcdb=sparse((-a.value./(a.value.^2+b.value.^2)).'(:));
            elseif sizeab(1)
              dcdb=sparse(numelc,numelb);
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(-a.value./(a.value.^2+b.value.^2)).'(:);
            else
              dcdb=sparse(numelc,numelb);
              dcdb(logical(repmat(speye(numelb),sizea(1),1)))=(-a.value./(a.value.^2+b.value.^2))(:);
            end
            dcda=spdiags((b.value./(a.value.^2+b.value.^2)).'(:),0,numelc,numelc);
          else
            if all(~sizeab)
              dcda=sparse((b.value./(a.value.^2+b.value.^2)).'(:));
            elseif sizeab(1)
              dcda=sparse(numelc,numela);
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=(b.value./(a.value.^2+b.value.^2)).'(:);
            else
              dcda=sparse(numelc,numela);
              dcda(logical(repmat(speye(numela),sizeb(1),1)))=(b.value./(a.value.^2+b.value.^2))(:);
            end
            dcdb=spdiags((-a.value./(a.value.^2+b.value.^2)).'(:),0,numelc,numelc);
          end
        end
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(atan2(a.value,b));
        sizea=size(a.value);
        sizeb=size(b);
        numela=numel(a.value);
        numelb=numel(b);
        numelc=numel(c.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcda=spdiags((b./(a.value.^2+b.^2)).'(:),0,numelc,numelc);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            dcda=spdiags((b./(a.value.^2+b.^2)).'(:),0,numelc,numelc);
          else
            if all(~sizeab)
              dcda=sparse((b./(a.value.^2+b.^2)).'(:));
            elseif sizeab(1)
              dcda=sparse(numelc,numela);
              dcda(logical(repelem(eye(numela),sizeb(2),1)))=(b./(a.value.^2+b.^2)).'(:);
            else
              dcda=sparse(numelc,numela);
              dcda(logical(repmat(speye(numela),sizeb(1),1)))=(b./(a.value.^2+b.^2))(:);
            end
          end
        end
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(atan2(a,b.value));
        sizea=size(a);
        sizeb=size(b.value);
        numela=numel(a);
        numelb=numel(b.value);
        numelc=numel(c.value);
        sizeab=sizea==sizeb;
        if all(sizeab)
          dcdb=spdiags((-a./(a.^2+b.value.^2)).'(:),0,numelc,numelc);
        else
          if any(sizea(~sizeab)>sizeb(~sizeab))
            if all(~sizeab)
              dcdb=sparse((-a./(a.^2+b.value.^2)).'(:));
            elseif sizeab(1)
              dcdb=sparse(numelc,numelb);
              dcdb(logical(repelem(eye(numelb),sizea(2),1)))=(-a./(a.^2+b.value.^2)).'(:);
            else
              dcdb=sparse(numelc,numelb);
              dcdb(logical(repmat(speye(numelb),sizea(1),1)))=(-a./(a.^2+b.value.^2))(:);
            end
          else
            dcdb=spdiags((-a./(a.^2+b.value.^2)).'(:),0,numelc,numelc);
          end
        end
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="atan2";
    endfunction

    function c=sind(a)
      numela=numel(a.value);
      c=ad(sind(a.value));
      dcda=spdiags(pi/180*cosd(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="sind";
    endfunction

    function c=cosd(a)
      numela=numel(a.value);
      c=ad(cosd(a.value));
      dcda=spdiags(-pi/180*sind(a.value).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cosd";
    endfunction

    function c=tand(a)
      numela=numel(a.value);
      c=ad(tand(a.value));
      dcda=spdiags(pi/180*(1+tand(a.value).^2).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="tand";
    endfunction

    function c=secd(a)
      numela=numel(a.value);
      c=ad(secd(a.value));
      dcda=spdiags((pi/180*tand(a.value).*secd(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="secd";
    endfunction

    function c=cscd(a)
      numela=numel(a.value);
      c=ad(cscd(a.value));
      dcda=spdiags((-pi/180*cotd(a.value).*cscd(a.value)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cscd";
    endfunction

    function c=cotd(a)
      numela=numel(a.value);
      c=ad(cotd(a.value));
      dcda=spdiags(pi/180*(-cotd(a.value).^2-1).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="cotd";
    endfunction

    function c=asind(a)
      numela=numel(a.value);
      c=ad(asind(a.value));
      dcda=spdiags((180/pi./(1-a.value.^2).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asind";
    endfunction

    function c=acosd(a)
      numela=numel(a.value);
      c=ad(acosd(a.value));
      dcda=spdiags((-180/pi./(1-a.value.^2).^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acosd";
    endfunction

    function c=atand(a)
      numela=numel(a.value);
      c=ad(atand(a.value));
      dcda=spdiags((180/pi./(1+a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="atand";
    endfunction

    function c=asecd(a)
      numela=numel(a.value);
      c=ad(asecd(a.value));
      dcda=spdiags((180/pi./(a.value.^2.*(1-1./a.value.^2).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="asecd";
    endfunction

    function c=acscd(a)
      numela=numel(a.value);
      c=ad(acscd(a.value));
      dcda=spdiags((-180/pi./(a.value.^2.*(1-1./a.value.^2).^.5)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acscd";
    endfunction

    function c=acotd(a)
      numela=numel(a.value);
      c=ad(acotd(a.value));
      dcda=spdiags((-180/pi./(a.value.^2+1)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="acotd";
    endfunction

    function c=atan2d(a,b)
      c=rad2deg(atan2(a,b));
      c.opfunc="atan2d";
    endfunction

    function c=erf(a)
      numela=numel(a.value);
      c=ad(erf(a.value));
      dcda=spdiags((2/pi^.5*exp(-a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="erf";
    endfunction

    function c=erfc(a)
      numela=numel(a.value);
      c=ad(erfc(a.value));
      dcda=spdiags((-2/pi^.5*exp(-a.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="erfc";
    endfunction

    function c=erfcx(a)
      numela=numel(a.value);
      c=ad(erfcx(a.value));
      dcda=spdiags((2*a.value.*exp(a.value.^2).*erfc(a.value)-2/pi^.5).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="erfcx";
    endfunction

    function c=erfinv(a)
      numela=numel(a.value);
      c=ad(erfinv(a.value));
      dcda=spdiags((pi^.5/2*exp(c.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="erfinv";
    endfunction

    function c=erfcinv(a)
      numela=numel(a.value);
      c=ad(erfcinv(a.value));
      dcda=spdiags((-pi^.5/2*exp(c.value.^2)).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="erfcinv";
    endfunction

    function c=ceil(a)
      numela=numel(a.value);
      c=ad(ceil(a.value));
      dcda=sparse(numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="ceil";
    endfunction

    function c=fix(a)
      numela=numel(a.value);
      c=ad(fix(a.value));
      dcda=sparse(numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="fix";
    endfunction

    function c=floor(a)
      numela=numel(a.value);
      c=ad(floor(a.value));
      dcda=sparse(numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="floor";
    endfunction

    function c=round(a)
      numela=numel(a.value);
      c=ad(round(a.value));
      dcda=sparse(numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="round";
    endfunction

    function c=abs(a)
      numela=numel(a.value);
      c=ad(abs(a.value));
      n=(real(a.value).^2+imag(a.value).^2).^.5;
      dcda=spdiags((real(a.value)./n).'(:),0,numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="abs";
    endfunction

    function c=arg(a)
      c=atan2(imag(a),real(a));
      c.opfunc="arg";
    endfunction

    function c=sign(a)
      c=a./abs(a);
      c.opfunc="sign";
    endfunction

    function c=real(a)
      numela=numel(a.value);
      c=ad(real(a.value));
      dcda=speye(numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="real";
    endfunction

    function c=imag(a)
      numela=numel(a.value);
      c=ad(imag(a.value));
      dcda=sparse(numela,numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="imag";
    endfunction

    function c=conj(a)
      numela=numel(a.value);
      c=ad(conj(a.value));
      dcda=speye(numela);
      c=insertdif1(c,a,dcda);
      c.opfunc="conj";
    endfunction

    function [W,IW]=max(X,Y,DIM)
      if nargin()==1
        sizeX=size(X);
        DIM=find(sizeX>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      if nargin()==2
        if isa(X,"ad")
          xv=X.value;
        else
          xv=X;
        end
        if isa(Y,"ad")
          yv=Y.value;
        else
          yv=Y;
        end
        sizex=size(X);
        sizey=size(Y);
        sizew=max(sizex,sizey);
        if all(sizex<sizew)
          xr=repmat(X,sizew);
        elseif sizex(1)<sizew(1)
          xr=repmat(X,sizew(1),1);
        elseif sizex(2)<sizew(2)
          xr=repmat(X,1,sizew(2));
        else
          xr=X;
        end
        if all(sizey<sizew)
          yr=repmat(Y,sizew);
        elseif sizey(1)<sizew(1)
          yr=repmat(Y,sizew(1),1);
        elseif sizey(2)<sizew(2)
          yr=repmat(Y,1,sizew(2));
        else
          yr=Y;
        end
        idx.type="()";
        w=max(xv,yv);
        W=ad(zeros(sizew));
        idx.subs={w==xv};
        if any(idx.subs{1}(:))
          W=subsasgn(W,idx,subsref(xr,idx));
        elseif isa(X,"ad")
          W=insertdif1(W,X,zeros(numel(W),numel(X)));
        end
        idx.subs={w==yv};
        if any(idx.subs{1}(:))
          W=subsasgn(W,idx,subsref(yr,idx));
        elseif isa(Y,"ad")
          W=insertdif1(W,Y,zeros(numel(W),numel(Y)));
        end
      else
        [~,IW]=max(X.value,[],DIM);
        idx.type="()";
        if DIM==1
          idx.subs={sub2ind(size(X.value),IW,reshape(1:size(X.value,(1:2)((1:2)~=DIM)),size(IW)))};
        else
          idx.subs={sub2ind(size(X.value),reshape(1:size(X.value,(1:2)((1:2)~=DIM)),size(IW)),IW)};
        end
        W=subsref(X,idx);
      end
      W.opfunc="max";
    endfunction

    function [W,IW]=min(X,Y,DIM)
      if nargin()==1
        sizeX=size(X);
        DIM=find(sizeX>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      if nargin()==2
        if isa(X,"ad")
          xv=X.value;
        else
          xv=X;
        end
        if isa(Y,"ad")
          yv=Y.value;
        else
          yv=Y;
        end
        sizex=size(X);
        sizey=size(Y);
        sizew=max(sizex,sizey);
        if all(sizex<sizew)
          xr=repmat(X,sizew);
        elseif sizex(1)<sizew(1)
          xr=repmat(X,sizew(1),1);
        elseif sizex(2)<sizew(2)
          xr=repmat(X,1,sizew(2));
        else
          xr=X;
        end
        if all(sizey<sizew)
          yr=repmat(Y,sizew);
        elseif sizey(1)<sizew(1)
          yr=repmat(Y,sizew(1),1);
        elseif sizey(2)<sizew(2)
          yr=repmat(Y,1,sizew(2));
        else
          yr=Y;
        end
        idx.type="()";
        w=min(xv,yv);
        W=ad(zeros(sizew));
        idx.subs={w==xv};
        if any(idx.subs{1}(:))
          W=subsasgn(W,idx,subsref(xr,idx));
        elseif isa(X,"ad")
          W=insertdif1(W,X,zeros(numel(W),numel(X)));
        end
        idx.subs={w==yv};
        if any(idx.subs{1}(:))
          W=subsasgn(W,idx,subsref(yr,idx));
        elseif isa(Y,"ad")
          W=insertdif1(W,Y,zeros(numel(W),numel(Y)));
        end
      else
        [~,IW]=min(X.value,[],DIM);
        idx.type="()";
        if DIM==1
          idx.subs={sub2ind(size(X.value),IW,reshape(1:size(X.value,(1:2)((1:2)~=DIM)),size(IW)))};
        else
          idx.subs={sub2ind(size(X.value),reshape(1:size(X.value,(1:2)((1:2)~=DIM)),size(IW)),IW)};
        end
        W=subsref(X,idx);
      end
      W.opfunc="min";
    endfunction

    function c=sumsq(a,DIM)
      if nargin<2||isempty(DIM)||length(DIM)~=1
        sizea=size(a.value);
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      if iscomplex(a.value)
        c=sum(a.*conj(a),DIM);
      else
        c=sum(a.*a,DIM);
      end
      c.opfunc="sumsq";
    endfunction

    function c=dot(a,b,DIM)
      if nargin<3||isempty(DIM)||length(DIM)~=1
        sizea=size(a);
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      if iscomplex(a)
        c=sum(conj(a).*b,DIM);
      else
        c=sum(a.*b,DIM);
      end
      c.opfunc="dot";
    endfunction

    function c=cross(a,b,DIM)
      if ~all(size(a)==size(b))
        error("cross: X and Y must have the same dimensions");
      end
      sizea=size(a);
      if nargin==3&&sizea(DIM)~=3
        error("cross: dimension DIM must have 3 elements");
      end
      if nargin<3||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea==3,1);
        if isempty(DIM)
          error("cross: must have at least one dimension with 3 elements");
        end
      end
      if DIM==1
        i1.type="()";
        i1.subs={1,":"};
        i2.type="()";
        i2.subs={2,":"};
        i3.type="()";
        i3.subs={3,":"};
        c=[subsref(a,i2).*subsref(b,i3)-subsref(a,i3).*subsref(b,i2);
        -(subsref(a,i1).*subsref(b,i3)-subsref(a,i3).*subsref(b,i1));
        subsref(a,i1).*subsref(b,i2)-subsref(a,i2).*subsref(b,i1)];
      elseif DIM==2
        i1.type="()";
        i1.subs={":",1};
        i2.type="()";
        i2.subs={":",2};
        i3.type="()";
        i3.subs={":",3};
        c=[subsref(a,i2).*subsref(b,i3)-subsref(a,i3).*subsref(b,i2) -(subsref(a,i1).*subsref(b,i3)-subsref(a,i3).*subsref(b,i1)) subsref(a,i1).*subsref(b,i2)-subsref(a,i2).*subsref(b,i1)];
      end
      c.opfunc="cross";
    endfunction

    function c=mean(a,varargin)
      DIM=[];
      opt="a";
      sizea=size(a.value);
      for i=1:length(varargin)
        if isa(varargin{i},"numeric")
          DIM=varargin{i};
        elseif strcmp(varargin{i},"a")||strcmp(varargin{i},"g")||strcmp(varargin{i},"h")
          opt=varargin{i};
        end
      end
      if isempty(DIM)
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      switch opt
        case "a"
          c=sum(a,DIM)/sizea(DIM);
        case "g"
          c=(-1).^(sum(a.value<0,DIM)./sizea(DIM)).*exp(sum(log(abs(a)),DIM)./sizea(DIM));
        case "h"
          c=sizea(DIM)./sum(a.^-1,DIM);
      end
      c.opfunc="mean";
    endfunction

    function c=median(a,DIM)
      sizea=size(a.value);
      if nargin<2||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      S=sort(a,DIM);
      idx.type="()";
      if rem(sizea(DIM),2)==0
        if DIM==1
          idx.subs={sizea(DIM)/2+(0:1),":"};
        elseif DIM==2
          idx.subs={":",sizea(DIM)/2+(0:1)};
        end
        c=sum(subsref(S,idx),DIM)/2;
      else
        if DIM==1
          idx.subs={ceil(sizea(DIM)/2),":"};
        elseif DIM==2
          idx.subs={":",ceil(sizea(DIM)/2)};
        end
        c=subsref(S,idx);
      end
      c.opfunc="median";
    endfunction

    function c=var(a,opt,DIM)
      sizea=size(a.value);
      if nargin<3||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      if nargin()<2||isempty(opt)||length(opt)~=1
        opt=0;
      end
      if sizea(DIM)==1
        opt=1;
      end
      switch opt
        case 0
          c=1/(sizea(DIM)-1)*sum((a-mean(a,DIM)).^2,DIM);
        case 1
          c=1/(sizea(DIM))*sum((a-mean(a,DIM)).^2,DIM);
      end
      c.opfunc="var";
    endfunction

    function c=std(a,varargin)
      c=sqrt(var(a,varargin{:}));
      c.opfunc="std";
    end

    %sums and products
    function c=sum(a,DIM)
      if nargin<2||isempty(DIM)||length(DIM)~=1
        sizea=size(a.value);
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      c=ad(sum(a.value,DIM));
      sizea=size(a);
      if DIM==1
        dcda=repmat(speye(sizea(2)),1,sizea(1));
      else
        dcda=sparse(repelem(eye(sizea(1)),1,sizea(2)));
      end
      c=insertdif1(c,a,dcda);
      c.opfunc="sum";
    endfunction

    function c=cumsum(a,DIM)
      sizea=size(a.value);
      if nargin<2||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      c=ad(cumsum(a.value,DIM));
      if DIM==1
        dcda=tril(repmat(speye(sizea(2)),sizea(1)));
      else
        dcda=sparse(tril(repelem(eye(sizea(1)),sizea(2),sizea(2))));
      end
      c=insertdif1(c,a,dcda);
      c.opfunc="cumsum";
    endfunction

    function c=prod(a,DIM)
      sizea=size(a.value);
      if nargin<2||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      c=ad(prod(a.value,DIM));
      dcda=sparse(numel(c.value),prod(sizea));
      numela=numel(a.value);
      if DIM==1
        vals=repmat(a.value,1,sizea(1));
        vals(logical(repelem(eye(sizea(1)),1,sizea(2))))=1;
        vals=prod(vals,1);
        dcda(logical(repmat(speye(sizea(2)),1,sizea(1))))=vals;
      else
        vals=repelem(a.value,sizea(2),1);
        vals(logical(repmat(eye(sizea(2)),sizea(1),1)))=1;
        vals=prod(vals,2);
        dcda(logical(repelem(eye(sizea(1)),1,sizea(2))))=vals;
      end
      c=insertdif1(c,a,dcda);
      c.opfunc="prod";
    endfunction

    function c=cumprod(a,DIM)
      sizea=size(a.value);
      if nargin<2||isempty(DIM)||length(DIM)~=1
        DIM=find(sizea>1,1);
        if isempty(DIM)
          DIM=1;
        end
      end
      numela=numel(a);
      c=ad(cumprod(a.value,DIM));
      dcda=sparse(numela,numela);
      if DIM==1
        vals=repmat(repelem(a.value,1,sizea(1)),1,sizea(1));
        vals(logical(repelem(eye(sizea(1)),1,sizea(1)*sizea(2))))=1;
        vals(logical(repmat(tril(ones(sizea(1)),-1),1,numela)))=1;
        vals=prod(vals,1);
        dcda(logical(repmat(speye(sizea(2)),sizea(1))))=vals;
        dcda=tril(dcda);
      else
        vals=repelem(a.value,sizea(2)^2,1);
        vals(logical(repmat(repelem(eye(sizea(2)),sizea(2),1),sizea(1),1)))=1;
        vals(logical(repmat(triu(ones(sizea(2)),1),numela,1)))=1;
        vals=prod(vals,2);
        dcda(logical(repelem(eye(sizea(1)),sizea(2),sizea(2))))=vals;
        dcda=tril(dcda);
      end
      c=insertdif1(c,a,dcda);
      c.opfunc="cumprod";
    endfunction

    %matrix functions
    function c=inv(a)
      c=ad(inv(a.value));
      dcda=kron(-c.value,c.value.');
      c=insertdif1(c,a,dcda);
      c.opfunc="inv";
    endfunction

    function c=pinv(a,tol)
      [U S V]=svd(a);
      diagS=diag(S);
      absdiagS=abs(diagS.value);
      if nargin()<2
        tol=max(size(a.value))*max(absdiagS)*eps;
      end
      idx.type="()";
      idx.subs={absdiagS<tol};
      diagS=subsasgn(diagS,idx,subsref(diagS,idx));
      idx.subs={absdiagS>=tol};
      diagS=subsasgn(diagS,idx,1./subsref(diagS,idx));
      S=diag(diagS);
      c=(U*S*V').';
      c.opfunc="pinv";
    endfunction

    function c=diag(a,varargin)
      c=ad(diag(a.value,varargin{:}));
      dcda=ad.difreshape(size(a.value),numel(a.value),numel(c.value),@(x)diag(x,varargin{:}));
      c=insertdif1(c,a,dcda);
      c.opfunc="diag";
    endfunction

    function c=trace(a)
      if (ndims(a.value)>2)
        error ("trace: only valid on 2-D cects");
      elseif (isempty(a.value))
        c=ad(0);
        c=insertdif1(c,a,sparse(numel(a.value)));
      elseif(isvector(a.value))
        c=ad(a.value(1));
        c=insertdif1(c,a,[1 sparse(1,numel(a.value)-1)]);
      else
        c=sum(diag(a));
      end
      c.opfunc="trace";
    endfunction

    function varargout=svd(a)
      nout=nargout();
      [U S V]=svd(a.value,"econ");
      m=size(a.value,1);
      n=size(a.value,2);
      diagS=diag(S);
      k=length(diagS);
      eyek=eye(k);
      F=1./(repelem(diagS.'.^2,k,1)-repelem(diagS.^2,1,k));
      F(sub2ind([k k],1:k,1:k))=0;
      if nout<=1
        c=ad(diag(S));
        dS=ad.difreshape([m n],m*n,k,@diag)*diag(speye(k).'(:))*kron(U',V.');
        c=insertdif1(c,a,dS);
        c.opfunc="svd";
        varargout={c};
      else
        dU=sparse(kron(U,eye(k))*(diag(F.'(:))*(kron(U',(V*S)')+ad.boxproduct(S*V',U')))+kron((eye(m)-U*U'),(V/S)'));
        dS=sparse(diag(eye(k).'(:))*kron(U',V.'));
        dV=sparse(kron(V,eye(k))*(diag(F.'(:))*(kron(S*U',V')+ad.boxproduct(V',(U*S)')))+ad.boxproduct((eye(n)-V*V'),(U/S)'));
        cdU=ad(U);
        cdU=insertdif1(cdU,a,dU);
        cdU.opfunc="svd";
        cdS=ad(S);
        cdS=insertdif1(cdS,a,dS);
        cdS.opfunc="svd";
        cdV=ad(V);
        cdV=insertdif1(cdV,a,dV);
        cdV.opfunc="svd";
        varargout={cdU,cdS,cdV};
      end
    endfunction

    function varargout=qr(a,b,opt)
      nout=nargout();
      m=size(a.value,1);
      n=size(a.value,2);
      if nout==3
        [q,r,p]=qr(a.value);
        ap=a*p;
      else
        [q,r]=qr(a.value);
      end
      sizea=size(a.value);
      sizeq=size(q);
      numela=numel(a.value);
      numelq=numel(q);
      ##da=dq*r+q*dr
      ##dq'*q+q'*dq=0
      eq=[ad.boxproduct(eye(sizeq),q')+kron(q',speye(sizeq)) sparse(numelq,numela)
          kron(speye(sizeq),r.') kron(q,speye(sizea(2)))];
      rhs=[sparse(numelq,numela);speye(numela)];

      qi=reshape(1:numelq,circshift(sizeq,1)).';
      qli=tril(qi,-1);
      qli(qli==0)=[];
      rli=tril(reshape(1:numela,circshift(sizea,1))',-1);
      rli(rli==0)=[];
      eq(qli,:)=[];
      rhs(qli,:)=[];
      eq(:,numelq+rli)=[];
      dqr=eq\rhs;
      dq=dqr(1:numelq,:);
######      dr=dqr(numelq+1:end,:);
      dr=sparse(numela,numela);
      rui=triu(reshape(1:numela,circshift(sizea,1)).').'(:);
      rui(rui==0)=[];
      dr(rui,:)=dqr(numelq+(1:length(rui)),:);

      if nout>1&&((nargin()==2&&b==0)||(nargin()==3&&opt==0))
        q(:,min(sizea)+1:end)=[];
        r(min(sizea)+1:end,:)=[];
      end
      QQ=ad(q);
      RR=ad(r);
      QQ.opfunc="qr";
      RR.opfunc="qr";
      if nout==3
        QQ=insertdif1(QQ,ap,dq);
        RR=insertdif1(RR,ap,dr);
      else
        QQ=insertdif1(QQ,a,dq);
        RR=insertdif1(RR,a,dr);
      end
      if nout==1
        varargout={RR};
      elseif nout==2
        if nargin()==1||(nargin()==2&&(b==0||isa(b,'char')))
          varargout={QQ,RR};
        else
          CC=QQ'*b;
          varargout={CC,RR};
        end
      elseif nout==3
        if ((nargin()==2&&(b==0||strcmp(b,'vector')))||(nargin()==3&&(opt==0||strcmp(opt,'vector'))))
          [~,p]=max(p,[],1);
        end
        varargout={QQ,RR,p};
      end
    endfunction

    function varargout=eig(a,b)
      nout=nargout();
      m=size(a,1);
      generlized=(nargin()==2&&isnumeric(b)&&all(size(a)==size(b)));
      if generlized
        if isa(a,"ad")&&isa(b,"ad")
          av=a.value;
          bv=b.value;
        elseif isa(a,"ad")
          av=a.value;
          bv=b;
        elseif isa(b,"ad")
          av=a;
          bv=b.value;
        end
      else
        av=a.value;
        bv=speye(m);
      end
      if nout<=2
        if generlized
          [v,l]=eig(av,bv);
        else
          [v,l]=eig(av);
        end
      else
        if generlized
          [v,l,w]=eig(av,bv);
        else
          [v,l,w]=eig(av);
        end
      end

      if isa(a,"ad")
        dl=sparse(diag(speye(m)(:))*kron(inv(bv*v),v.'));
        if nout>1
          eq=sparse(kron(bv,l))-kron(av,speye(m));
          rhs=((kron(speye(m),v')-kron(bv*v,speye(m))*dl));
          if isreal(v)&&generlized
            nor=sparse(m,m^2);
            nor(sub2ind([m,m^2],(1:m)(:),reshape(1:m^2,m,m)'(abs(v)==max(abs(v),[],1))))=1;
          else
            nor=repmat(speye(m),1,m)*spdiags(v'(:),0,m^2,m^2);
          end
          eq=[eq;nor];
          rhs=[rhs;sparse(m,m^2)];
          dv=(eq\rhs);
          if nout==3
            bp=sparse(ad.boxproduct(eye(m),eye(m)));
            eq=sparse(kron(bv',l))-kron(av',speye(m));
            rhs=((kron(speye(m),w')-kron(bv'*w,speye(m))*dl*bp));
            if generlized
              nor=sparse(m,m^2);
              nor(sub2ind([m,m^2],(1:m)(:),reshape(1:m^2,m,m)'(abs(w)==max(abs(w),[],1))))=1;
            else
              nor=repmat(speye(m),1,m)*spdiags(w'(:),0,m^2,m^2);
            end
            eq=[eq;nor];
            rhs=[rhs;sparse(m,m^2)];
            dw=(eq\rhs)*bp;
          end
        end
      end
      if generlized&&isa(b,"ad")
        eq=[kron(av,speye(m))-sparse(kron(bv,l)) -kron(bv*v,speye(m))*diag(eye(m)(:))];
        rhs=[((kron(speye(m),(v*l)')))];
        if isreal(v)
          nor=sparse(m,m^2);
          nor(sub2ind([m,m^2],(1:m)(:),reshape(1:m^2,m,m)'(abs(v)==max(abs(v),[],1))))=1;
        else
          nor=repmat(speye(m),1,m)*spdiags(v'(:),0,m^2,m^2);
        end
        eq=[eq;nor sparse(m,m^2)];
        rhs=[rhs;sparse(m,m^2)];
        dlvdb=(eq\rhs);
        dvdb=dlvdb(1:m^2,:);
        dldb=dlvdb(m^2+(1:m^2),:);
        if generlized&&isa(b,"ad")
          bp=sparse(ad.boxproduct(eye(m),eye(m)));
          eq=kron(av',speye(m))-sparse(kron(bv',l));
          rhs=((kron(speye(m),(w*l)')+kron(bv'*w,speye(m))*dldb*bp));
          if isreal(w)
            nor=sparse(m,m^2);
            nor(sub2ind([m,m^2],(1:m)(:),reshape(1:m^2,m,m)'(abs(w)==max(abs(w),[],1))))=1;
          else
            nor=repmat(speye(m),1,m)*spdiags(w'(:),0,m^2,m^2);
          end
          eq=[eq;nor];
          rhs=[rhs;zeros(m,m^2)];
          dwdb=(eq\rhs)*bp;
        end
      end

      lambda=ad(l);
      if generlized&&isa(b,"ad")&&isa(a,"ad")
        lambda=insertdif2(lambda,a,b,dl,dldb);
      elseif generlized&&isa(b,"ad")
        lambda=insertdif1(lambda,b,dldb);
      else
        lambda=insertdif1(lambda,a,dl);
      end
      lambda.opfunc="eig";
      if nout>1
        vec=ad(v);
        if generlized&&isa(b,"ad")&&isa(a,"ad")
          vec=insertdif2(vec,a,b,dv,dvdb);
        elseif generlized&&isa(b,"ad")
          vec=insertdif1(vec,b,dvdb);
        else
          vec=insertdif1(vec,a,dv);
        end
        vec.opfunc="eig";
      end
      if nout<=1
        varargout={diag(lambda)};
      elseif nout==2
        varargout={vec,lambda};
      elseif nout==3
        wvec=ad(w);
        if generlized&&isa(b,"ad")&&isa(a,"ad")
          wvec=insertdif2(wvec,a,b,dw,dwdb);
        elseif generlized&&isa(b,"ad")
          wvec=insertdif1(wvec,b,dwdb);
        else
          wvec=insertdif1(wvec,a,dw);
        end
        wvec.opfunc="eig";
        varargout={vec,lambda,wvec};
      end
    endfunction

    function varargout=eigs(a,varargin)
      nout=nargout();
      m=size(a,1);
      generlized=length(varargin)>0&&isnumeric(varargin{1})&&all(size(varargin{1})==size(a));
      if generlized
        if isa(varargin{1},"ad")
          bv=varargin{1}.value;
        else
          bv=varargin{1};
        end
        if isa(a,"ad")
          av=a.value;
        else
          av=a;
        end
        [v,l,f]=eigs(av,bv,varargin{2:end});
      else
        av=a.value;
        bv=eye(m);
        [v,l,f]=eigs(av,varargin{:});
      end

      K=size(l,1);
      if isa(a,"ad")
        eq=[sparse(kron(bv,l))-kron(av,speye(K)) kron(bv*v,speye(K))*diag(eye(K)(:))];
        rhs=((kron(speye(m),v')));
        if isreal(v)&&generlized
          nor=sparse(K,K*m);
          nor(sub2ind([K,K*m],(1:K)(:),reshape(1:K*m,K,m)'(abs(v)==max(abs(v),[],1))))=1;
        else
          nor=repmat(speye(K),1,m)*spdiags(v'(:),0,K*m,K*m);
        end
        eq=[eq;nor sparse(K,K^2)];
        rhs=[rhs;sparse(K,m^2)];
        dvl=(eq\rhs);
        dv=dvl(1:m*K,:);
        dl=dvl(m*K+(1:K^2),:);
      end
      if generlized&&isa(varargin{1},"ad")
        eq=[kron(av,speye(K))-sparse(kron(bv,l)) -kron(bv*v,speye(K))*diag(eye(K)(:))];
        rhs=[((kron(speye(m),(v*l)')))];
        if isreal(v)&&generlized
          nor=sparse(K,K*m);
          nor(sub2ind([K,K*m],(1:K)(:),reshape(1:K*m,K,m)'(abs(v)==max(abs(v),[],1))))=1;
        else
          nor=repmat(speye(K),1,m)*spdiags(v'(:),0,K*m,K*m);
        end
        eq=[eq;nor sparse(K,K^2)];
        rhs=[rhs;sparse(K,m^2)];
        dlvdb=(eq\rhs);
        dvdb=dlvdb(1:m*K,:);
        dldb=dlvdb(m*K+(1:K^2),:);
      end

      lambda=ad(l);
      if generlized&&isa(varargin{1},"ad")&&isa(a,"ad")
        lambda=insertdif2(lambda,a,varargin{1},dl,dldb);
      elseif generlized&&isa(varargin{1},"ad")
        lambda=insertdif1(lambda,varargin{1},dldb);
      else
        lambda=insertdif1(lambda,a,dl);
      end
      lambda.opfunc="eigs";
      if nout>1
        vec=ad(v);
        if generlized&&isa(varargin{1},"ad")&&isa(a,"ad")
          vec=insertdif2(vec,a,varargin{1},dv,dvdb);
        elseif generlized&&isa(varargin{1},"ad")
          vec=insertdif1(vec,varargin{1},dvdb);
        else
          vec=insertdif1(vec,a,dv);
        end
        vec.opfunc="eigs";
      end
      if nout==1
        varargout={diag(lambda)};
      elseif nout==2
        varargout={vec,lambda};
      elseif nout==3
        varargout={vec,lambda,f};
      end
    endfunction

    function c=kron(a,b)
      if isa(a,"ad")&&isa(b,"ad")
        c=ad(kron(a.value,b.value));
        sizea=size(a.value);
        sizeb=size(b.value);
        sizec=size(c.value);
        numela=numel(a.value);
        numelb=numel(b.value);
        numelc=numel(c.value);
        dcda=sparse(numelc,numela);
        dcda(logical(ad.difreshape(sizea,numela,numelc,@(x)kron(x,ones(sizeb)))))=repmat(b.value.'(:),numela,1);
        dcdb=sparse(numelc,numelb);
        dcdb(logical(ad.difreshape(sizeb,numelb,numelc,@(x)kron(ones(sizea),x))))=repmat(a.value.'(:),numelb,1);
        c=insertdif2(c,a,b,dcda,dcdb);
      elseif isa(a,"ad")
        c=ad(kron(a.value,b));
        sizea=size(a.value);
        sizeb=size(b);
        sizec=size(c.value);
        numela=numel(a.value);
        numelc=numel(c.value);
        dcda=sparse(numelc,numela);
        dcda(logical(ad.difreshape(sizea,numela,numelc,@(x)kron(x,ones(sizeb)))))=repmat(b.'(:),numela,1);
        c=insertdif1(c,a,dcda);
      elseif isa(b,"ad")
        c=ad(kron(a,b.value));
        sizea=size(a);
        sizeb=size(b.value);
        sizec=size(c.value);
        numelb=numel(b.value);
        numelc=numel(c.value);
        dcdb=sparse(numelc,numelb);
        dcdb(logical(ad.difreshape(sizeb,numelb,numelc,@(x)kron(ones(sizea),x))))=repmat(a.'(:),numelb,1);
        c=insertdif1(c,b,dcdb);
      end
      c.opfunc="kron";
    endfunction

    function c=sqrtm(a)
      n=size(a.value,1);
      c=ad(sqrtm(a.value));
      dcda=sparse(inv(kron(eye(n),c.value.')+kron(c.value,eye(n))));
      c=insertdif1(c,a,dcda);
      c.opfunc="sqrtm";
    endfunction

    function c=det(a)
      n=size(a.value,1);
      c=ad(det(a.value));
      dcda=sparse((inv(a.value)*c.value)(:).');
##      da=difsizederiv(1,n^2,[n n],@(d)c.value*trace(a.value\d));
      c=insertdif1(c,a,dcda);
      c.opfunc="det";
    endfunction

    function [c,P]=chol(a,opt)
      if nargin()==1
        opt="upper";
      end
      nout=nargout();
      if nout==2
        [U P]=chol(a.value,opt);
        if P~=0
          c=U;
          warning('Cholesky factorization failed, no derivative returned.');
          return
        end
      else
        U=chol(a.value,opt);
      end
      m=size(a.value,1);
      if strcmp(opt,"upper")
        eq=ad.boxproduct(eye(m),U.')+kron(U.',eye(m));
        ueq=triu(reshape(1:m^2,m,[])')(:);
      else
        eq=kron(eye(m),conj(U))+ad.boxproduct(U,eye(m));
        ueq=tril(reshape(1:m^2,m,[])')(:);
      end
      ueq(ueq==0)=[];
      dcda=sparse(m^2,m^2);
      dcda(ueq,ueq)=eq(ueq,ueq)\(eye(m^2)(ueq,ueq));
      c=ad(U);
      c=insertdif1(c,a,dcda);
      c.opfunc="chol";
    endfunction

    function [L,U,P]=lu(a)
      [l u P]=lu(a.value);
      sizea=size(a.value);
      numela=numel(a.value);
      dl=zeros(numel(l),numela);
      du=zeros(numel(u),numela);
      i=1;
      j=1;
      ids=1:numela;
      ids=reshape(ids,sizea(2),[]).';
      trilids=tril(ids,-1)(:);
      trilids(trilids==0)=[];
      triuids=triu(ids);
      triuids(triuids==0)=[];

      lc=kron(speye(sizea(1)),[u;sparse(max(sizea(2)-sizea(1),0),sizea(2))].');
      uc=kron([l sparse(sizea(1),max(sizea(1)-sizea(2),0))],speye(sizea(2)));
      lc(triuids,triuids)=0;
      uc(trilids,trilids)=0;
      sol=(lc+uc)\speye(numela);
      dl=sparse(numela,numela);
      dl(trilids,:)=sol(trilids,:);
      du=sparse(numela,numela);
      du(triuids,:)=sol(triuids,:);
      if nargout()<=1
        if sizea(2)>sizea(1)
          l=[l,zeros(sizea(1),sizea(2)-sizea(1))];
        elseif m<n
          u=[u;zeros(sizea(1)-sizea(2),sizea(2))];
        end
      else
        dl(ids(:,(1:sizea(2))>sizea(1)),:)=[];
        du(ids((1:sizea(1))>sizea(2),:),:)=[];
      end

      L=ad(l);
      L=insertdif1(L,P*a,dl);
      L.opfunc="lu";
      U=ad(u);
      U=insertdif1(U,P*a,du);
      U.opfunc="lu";
      if nargout()==2
        L=P\L;
      elseif nargout()<=1
        L=L+U-eye(sizea);
      end
    endfunction

    function X=sylvester(A,B,C)
      Xsize=size(C.value);
      s.type="()";
      s.subs={":"};
      X=reshape((kron(A,eye(Xsize(2)))+kron(eye(Xsize(1)),B'))\subsref(C',s),circshift(Xsize,1))';
      X.opfunc="sylvester";
    endfunction

    function varargout=schur(a)
      [u s]=schur(a.value);
      n=size(a.value,1);
      ids=reshape(1:n^2,n,n)';
      trilids=tril(ids,-1)(:);
      trilids(trilids==0)=[];
      zeroids=ids(s==0);
      if isreal(s)
        cid=setdiff(trilids,zeroids);
        lencid=length(cid);
        eqc=sparse(lencid,2*n^2);
        eqc(sub2ind([lencid,2*n^2],(1:lencid)(:),(n^2+cid+1)(:)))=1;
        eqc(sub2ind([lencid,2*n^2],(1:lencid)(:),(n^2+cid-n)(:)))=-1;
      else
        cid=2;
        lencid=1;
        eqc=sparse(lencid,2*n^2);
        eqc(sub2ind([lencid,2*n^2],(1:lencid)(:),(n^2+cid+1)(:)))=1;
        eqc(sub2ind([lencid,2*n^2],(1:lencid)(:),(n^2+cid-n)(:)))=-1;
      end


      eq=[ad.boxproduct(eye(n),u.')+kron(u',speye(n)) sparse(n^2,n^2);
            kron(u',s.')+kron(-u'*a.value,speye(n)) kron(speye(n),speye(n));
            eqc];
      rhs=[sparse(n^2,n^2);kron(u',u.');sparse(lencid,n^2)];
      eq(n^2+zeroids,n^2+zeroids)=0;
      eq(trilids,:)=[];
      rhs(trilids,:)=[];
      duds=eq\rhs;
      du=duds(1:n^2,:);
      ds=duds(n^2+1:2*n^2,:);

      S=ad(s);
      S=insertdif1(S,a,ds);
      S.opfunc="schur";
      if nargout()==1
        varargout={S};
      else
        U=ad(u);
        U=insertdif1(U,a,du);
        U.opfunc="schur";
        varargout={U,S};
      end
    endfunction

    function [U S]=rsf2csf(U,S)
      N=size(S,1);
      for m=N:-1:2
        submm1.type="()";
        submm1.subs={m,m-1};
        subm1m1.type="()";
        subm1m1.subs={m-1,m-1};
        submm.type="()";
        submm.subs={m,m};
        if abs(subsref(S,submm1))>eps()*(abs(subsref(S,subm1m1))+abs(subsref(S,submm)))
          subm1m.type="()";
          subm1m.subs={[m-1:m],[m-1:m]};
          mu=eig(subsref(S,subm1m))-subsref(S,submm);
          sub1.type="()";
          sub1.subs={1};
          r=norm([subsref(mu,sub1) subsref(S,submm1)]);
          c=subsref(mu,sub1)/r;
          s=subsref(S,submm1)/r;
          G=[c' s;-s c];
          Gc=G';
          j=[m-1 N];
          subkj.type="()";
          subkj.subs={[m-1:m],[m-1:N]};
          S=subsasgn(S,subkj,G*subsref(S,subkj));
          subik.type="()";
          subik.subs={[1:m],[m-1:m]};
          S=subsasgn(S,subik,subsref(S,subik)*Gc);
          subik.type="()";
          subik.subs={[1:N],[m-1:m]};
          U=subsasgn(U,subik,subsref(U,subik)*Gc);
        end
        S=subsasgn(S,submm1,0);
      end
      U.opfunc="rsf2csf";
      S.opfunc="rsf2csf";
    endfunction

    function varargout=balance(a,varargin)
      nout=nargout();
      if length(varargin)<1||isa(varargin{1},"char")
        [d,p,aa]=balance(a.value,varargin{:});
        dd=eye(size(a.value))(:,p)*diag(d);
        A=dd\a*dd;
        A.opfunc="balance";
        if nout==1
          varargout={A};
        elseif nout==2
          varargout={dd,A};
        elseif nout==3
          varargout={d,p,A};
        end
      else
        if isa(a,"ad")&&isa(varargin{1},"ad")
          [cc,dd,aa,aa]=balance(a.value,varargin{1}.value,varargin{2:end});
          A=cc*a*dd;
          A.opfunc="balance";
          B=cc*varargin{1}*dd;
          B.opfunc="balance";
          varargout={cc,dd,A,B};
        elseif isa(a,"ad")
          [cc,dd,aa,aa]=balance(a.value,varargin{:});
          A=cc*a*dd;
          A.opfunc="balance";
          B=cc*varargin{1}*dd;
          varargout={cc,dd,A,B};
        elseif isa(varargin{1},"ad")
          [cc,dd,aa,aa]=balance(a,varargin{1}.value,varargin{2:end});
          A=cc*a*dd;
          B=cc*varargin{1}*dd;
          B.opfunc="balance";
          varargout={cc,dd,A,B};
        end
      end
    endfunction

    function v=norm(a,p,opt)
      if nargin()==1
        p=2;
      end
      if all(size(a)>1)&&(nargin()<3||~(strcmp(opt,"rows")||strcmp(opt,"cols")||strcmp(opt,"columns")))
        if p==1
          sum1=sum(abs(a),1);
          maxi=find(sum1.value==max(sum1.value),1);
          s.type="()";
          s.subs={maxi};
          v=subsref(sum1,s);
        elseif p==2
          sigv=svd(a);
          maxi=find(sigv.value==max(sigv.value),1);
          s.type="()";
          s.subs={maxi};
          v=subsref(sigv,s);
        elseif p==inf || strcmpi(p,"inf")
          sum2=sum(abs(a),2);
          maxi=find(sum2.value==max(sum2.value),1);
          s.type="()";
          s.subs={maxi};
          v=subsref(sum2,s);
        elseif strcmpi(p,"fro")
          v=sqrt(sum(diag(a'*a)));
        elseif p>1
          if isa(p,"ad")
            if isa(a,"ad")
              v=ad(norm(a.value,p.value));
              dvda=sparse(ad.fd(@(a)norm(a,p.value),a.value));
              dvdp=sparse(ad.fd(@(p)norm(a.value,p),p.value));
              v=insertdif2(v,a,p,dvda,dvdp);
            else
              v=ad(norm(a,p.value));
              dvdp=sparse(ad.fd(@(p)norm(a,p),p.value));
              v=insertdif1(v,p,dvdp);
            end
          else
            v=ad(norm(a.value,p));
            dvda=sparse(ad.fd(@(a)norm(a,p),a.value));
            v=insertdif1(v,a,dvda);
          end
          warning("norm: Finite difference used to find matrix p norm.");
        end
      else
        if nargin()==3
          if strcmp(opt,"rows")
            sumdir=2;
          elseif strcmp(opt,"cols")||strcmp(opt,"columns")
            sumdir=1;
          end
          if p==inf || strcmpi(p,"inf")
            [~,maxi]=max(abs(a.value),[],sumdir);
            s.type="()";
            if sumdir==1
              s.subs={sub2ind(size(a.value),maxi,reshape(1:size(a.value,(1:2)((1:2)~=sumdir)),size(maxi)))};
            else
              s.subs={sub2ind(size(a.value),reshape(1:size(a.value,(1:2)((1:2)~=sumdir)),size(maxi)),maxi)};
            end
            v=abs(subsref(a,s));
          elseif p==-inf
            [~,mini]=min(abs(a.value),[],sumdir)
            s.type="()";
            if sumdir==1
              s.subs={sub2ind(size(a.value),mini,reshape(1:size(a.value,(1:2)((1:2)~=sumdir)),size(mini)))};
            else
              s.subs={sub2ind(size(a.value),reshape(1:size(a.value,(1:2)((1:2)~=sumdir)),size(mini)),mini)};
            end
            v=abs(subsref(a,s));
          elseif strcmpi(p,"fro")
            v=sum(abs(a).^2,sumdir).^.5;
          elseif p==0
            v=ad(norm(a.value,0,opt));
            dvda=sparse(numel(v.value),numel(a.value));
            v=insertdif1(v,a,dvda);
          else
            v=sum(abs(a).^p,sumdir).^(1/p);
          end
        else
          if p==inf || strcmpi(p,"inf")
            maxi=find(abs(a.value)==max(abs(a.value)),1);
            s.type="()";
            s.subs={maxi};
            v=abs(subsref(a,s));
          elseif p==-inf
            maxi=find(abs(a.value)==min(abs(a.value)),1);
            s.type="()";
            s.subs={maxi};
            v=abs(subsref(a,s));
          elseif strcmpi(p,"fro")
            v=sum(abs(a).^2)^.5;
          elseif p==0
            v=ad(norm(a.value,0));
            dvda=sparse(1,numel(a.value));
            v=insertdif1(v,a,dvda);
          else
            v=sum(abs(a).^p)^(1/p);
          end
        end
      end
      v.opfunc="norm";
    endfunction

  endmethods
endclassdef
## Generates id for objects
function n=numtimescalled()
  persistent m
  if isempty(m)
    m=uint64(1);
  else
    m++;
  end
  n=sprintf('%d',m);
endfunction
## Not used but can layout derivative correctly
function d=difsizederiv(nfoutput,nvars,sizevars,derivf)
  d=zeros(nfoutput,nvars);
  i=1;
  j=1;
  for k=1:nvars
    d1=zeros(sizevars);
    d1(i,j)=1;
    d1=derivf(d1);
    d(:,k)=d1.'(:);
    j++;
    if j>sizevars(2)
      i++;
      j=1;
    end
  end
endfunction

%!test adtest
