%This script creates the output data structure from the raw data.
%The output contains the data for the 500 firms with the most limit
%violations
%For each firm, it contains the longest data segment with at most 5 tradign
%days of data gap
%for each firm output has the following fields
%    output(ifirms).PriceDates: Dates of observation in matlab date frmat
%    output(ifirms).Price: Closing Price
%    output(ifirms).MarketCap: Market Cap in millions RMB 
%    output(ifirms).timestamp: timestamp in trading days. This is important
%    as day count convention in model is defined in trading days not
%    calendar days, and there are data gaps in the data.

clear all
 
addpath '.\CensoredGARCH';
 
addpath '.\SMC';
 
randn('state',0);
rand('state',0);
 
addpath '.\dailyprice_mktcap';
 
%%%%%%%%%%%
%load data%
%%%%%%%%%%%

%sort firms by total number of limit violations in descending order

 [NUM,TXT,RAW]=xlsread('Mapping_Table.xlsx','sheet1');

 firmnames=RAW(2:end,1);
 
 sourcefile=RAW(2:end,6);
 
 ST_inout=RAW(2:end,8:13);
 
 ST_mark=RAW(2:end,7);
 
 NumViolations=sum(NUM(:,2:3),2);
 
 [NumViolations,j]=sort(NumViolations,'descend');

 firmnames=firmnames(j);
 
 sourcefile=sourcefile(j);
 
 ST_inout=ST_inout(j,:);
 
 ST_mark=ST_mark(j,:);
 
 %trading days
 [NUM,TXT,RAW]=xlsread('Mapping_Table.xlsx','TradingDate');
  
 Trading_Date=datenum(TXT(2:end),'mm/dd/yyyy');
 
%%%%%%%%%%%%%%
%extract data%
%%%%%%%%%%%%%%

%number of firms that are saved to output file
Nfirms=500;

%number of trading days with no data that is allowed in one data sample
allowedgap=5;

clear viol_5 viol_10 N;

for ifirms=1:Nfirms
    
    
    ifirms
    %load data
        
    [NUM,TXT,RAW]=xlsread([sourcefile{ifirms} '.xlsx'],'dailyprice');
    
    
     [NUM2,TXT2,RAW2]=xlsread([sourcefile{ifirms} '.xlsx'],'mktcap');
    
    Ncol=find(strcmp(RAW(1,:),firmnames{ifirms}));
    
    %price
    a=RAW(2:end,Ncol);
    
    j=find(strcmp(a,'#N/A N/A'));
    
    a(j)=mat2cell(nan(length(j),1),ones(length(j),1));
    
    Price=cell2mat(a);
        
    %market cap
    a=RAW2(2:end,Ncol);
    
    j=find(strcmp(a,'#N/A N/A'));
    
    a(j)=mat2cell(nan(length(j),1),ones(length(j),1));
    
    MarketCap=cell2mat(a);
           
    T=length(Price);
    
    PriceDates=RAW(2:end,1);
    
    IsDate=ones(size(PriceDates));
    
    for i=1:T
       
        if isnan(PriceDates{i})
            
            IsDate(i)=0;
        end
        
    end
    
    keep=find(IsDate);
        
    PriceDates=datenum(PriceDates(keep),'mm/dd/yyyy');
    
    Price=Price(keep);    
    MarketCap=MarketCap(keep);
    
    %only keep trading dates
    [ph,keep]=intersect(PriceDates,Trading_Date);
    
    PriceDates=PriceDates(keep);
    Price=Price(keep);
    MarketCap=MarketCap(keep);
        
    %compute for each time the number of uninterrupted periods before
    T=size(Price,1);
        
    PrKeep=[]; PrDateKeep=[]; MarketCapKeep=[];

    obsstamp=[];
    
    segmentid=[];
    
    gap=0; 
    
    segmentcount=1;
    
    for t=1:T
       
        if isfinite(Price(t))
           
            PrKeep=[PrKeep; Price(t)];
            PrDateKeep=[PrDateKeep; PriceDates(t)];
            MarketCapKeep=[MarketCapKeep; MarketCap(t)];
            
            segmentid=[segmentid; segmentcount];
            
            obsstamp=[obsstamp; t];
            
            if gap>allowedgap
                segmentcount=segmentcount+1;
            end
                        
            gap=0;  
            
        else
            
            gap=gap+1;
            
        end
        
    end    
    
    %keep the largest data segment (within each data segment the maximum number of
    %trading days with no data is allowedgap)
    ta=tabulate(segmentid);
        
    [N(ifirms),jmax]=max(ta(:,2));
    
    keep=find(segmentid==ta(jmax,1));
        
    output(ifirms).PriceDates=PrDateKeep(keep);
    
    output(ifirms).Price=PrKeep(keep);
        
    output(ifirms).MarketCap=MarketCapKeep(keep);
    
    output(ifirms).timestamp=obsstamp(keep);
            
    ret=output(ifirms).Price(2:end)./output(ifirms).Price(1:end-1)-1;
    
    viol_5(ifirms)=mean(abs(ret)>.05-.001 & abs(ret)<.05+.001);
    
    viol_10(ifirms)=mean(abs(ret)>.1-.001 & abs(ret)<.1+.001);
        
end

save ChineseData output viol_5 viol_10 N allowedgap firmnames sourcefile ST_inout ST_mark;
