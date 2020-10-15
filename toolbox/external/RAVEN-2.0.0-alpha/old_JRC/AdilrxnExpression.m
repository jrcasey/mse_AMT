function results=AdilrxnExpression(model,file)
% The aim of this function is to decide the reaction in the Model is
% upregulated, down regulated or unchanged reactions from the expression file
% model is GEM
% file is the fileName containing gene Pvalue Foldchange

[Ensemble PValue fold ]=textread(file,'%s%f%f');
for i=1:numel(model.rxns) 
    g=find(full(model.rxnGeneMat(i,:)));
    genes=model.genes(g);
    J=find(ismember(Ensemble,genes));
    results.rxns(i)=model.rxns(i);
    results.numGene(i)=numel(J);
    pv=PValue(J);
    f=fold(J);
    I=(pv<0.05);
    J_new = J(I);
    f1=f;
    if numel(I)==0
        results.numUn(i)=numel(pv);
        results.numPos(i)=0;
        results.numNeg(i)=0;
        results.express(i)=0;
    else
        f=f(I);
        results.numPos(i)=numel(find(f>=0));
        results.numNeg(i)=numel(find(f<0));
        results.numUn(i)=numel(pv)-numel(f);
        if results.numPos(i) >=1 && results.numNeg(i) ==0
           results.express(i)=1;
        elseif results.numNeg(i) >=1 && results.numPos(i) ==0
           results.express(i)=-1;
        else
          results.express(i) =0;
        end      
        s=sprintf('%s(p=%f,f=%f);',cell2mat(Ensemble(J_new(1))),pv(I(1)),f1(I(1)));
        for k=2:numel(I)
           s=sprintf('%s%s(p=%f,f=%f);',s,cell2mat(Ensemble(J_new(k))),pv(I(k)),f1(I(k)));
        end
           results.report(i)={s};
    end
end
    
