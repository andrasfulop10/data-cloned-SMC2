
clear all;

Bounds=[0.05 0.1];

load ChineseData;

violations=[];

for ifirms=1:100
    
    Ret=output(ifirms).Price(2:end)./output(ifirms).Price(1:end-1)-1;
        
    violations(ifirms).PriceDates=output(ifirms).PriceDates(2:end);
    
    violations(ifirms).index=[];
        
    %round return to bounds
    for iBound=1:length(Bounds)
   
        index=find(abs(Ret)>Bounds(iBound)-.001 & abs(Ret)<Bounds(iBound)+.001);
    
        Ret(index)=sign(Ret(index)).*Bounds(iBound);
    
        violations(ifirms).index=[violations(ifirms).index; index];
        
    end
       
    
    violations(ifirms).Ret=Ret;
    
    violations(ifirms).T=length(Ret);
    
    violations(ifirms).vr=length(violations(ifirms).index)/violations(ifirms).T;
        
    violations(ifirms).firstv=min(violations(ifirms).index)/violations(ifirms).T;
    
end

load ChineseRes_GARCH;

res_garch=res;

load ChineseRes_CensoredGARCH;

res_cgarch=res;

load ChineseData;

[ph,j]=sort([violations.firstv],'descend');

for i=1:1;

    figure(i);
    
subplot(2,1,1);

plot(violations(j(i)).PriceDates,violations(j(i)).Ret,'-');

hold on;

plot(violations(j(i)).PriceDates(violations(j(i)).index),violations(j(i)).Ret(violations(j(i)).index),'ro')

hold off;

datetick;

title(['Stock Returns of Guanghui Logistics']);

subplot(2,1,2);

plot(violations(j(i)).PriceDates,[res_garch(j(i)).FiltMeans' res_cgarch(j(i)).FiltMeans.h']);

title('Variance Series (Censored vs Uncensored)');

legend('NGARCH(1,1)','Censored NGARCH(1,1)','Location','NorthWest');

datetick;

end

print -dpdf Figure_CensoredGARCH

for i=1:100
    
   mh_cg(i)=mean(res_cgarch(i).FiltMeans.h);
    
   vh_cg(i)=var(res_cgarch(i).FiltMeans.h);
   
   
   mh_g(i)=mean(res_garch(i).FiltMeans);
      
   vh_g(i)=var(res_garch(i).FiltMeans);
         
  MarketCap(i)=nanmean(output(i).MarketCap);
  
end

[ph,j]=sort([violations.vr],'descend');

i1=1:33;

[mean(MarketCap(j(i1)))/1000  mean([violations(j(i1)).vr]) mean(mh_cg(j(i1)))/mean(mh_g(j(i1))) mean(vh_cg(j(i1)))/mean(vh_g(j(i1)))]

i1=34:67;

[mean(MarketCap(j(i1)))/1000 mean([violations(j(i1)).vr]) mean(mh_cg(j(i1)))/mean(mh_g(j(i1))) mean(vh_cg(j(i1)))/mean(vh_g(j(i1)))]

i1=68:100;

[mean(MarketCap(j(i1)))/1000 mean([violations(j(i1)).vr]) mean(mh_cg(j(i1)))/mean(mh_g(j(i1))) mean(vh_cg(j(i1)))/mean(vh_g(j(i1)))]



