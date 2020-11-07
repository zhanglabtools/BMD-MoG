function [model,k,R] = Ktuning(model,k,R)
Sigma=model.Sigma;
weight=model.weight;
a=size(Sigma,2);
TOL=1e-5;

for i=1:a-1
    for j=i+1:a
        diff=abs(Sigma(i)-Sigma(j))/(Sigma(i)+Sigma(j));
        if diff<TOL
            Sigma(i)=(Sigma(i)+Sigma(j))/2;
            Sigma(j)=0;
            weight(i)=weight(i)+weight(j);
            weight(j)=0;
            R(:,i)=R(:,i)+R(:,j);
            R(:,j)=0;
            k=k-1;
        end
    end
end
Sigma(Sigma==0)=[];
weight(weight==0)=[];
R(:,sum(R,1)==0)=[];
model.mu=zeros(1,k);
model.Sigma=Sigma;
model.weight=weight;