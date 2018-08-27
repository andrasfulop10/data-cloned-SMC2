  function [logl,mtx_bin]=Resample(lnw,gpu_switch,R)
  %compute loglikelihood from loglikelihoods columnwise (from lnw)
  %also return R *Nparam random indices that directly index into lnw using
  %stratified resampling
  
  [Nparticles,Nparam]=size(lnw);
        
  lnw(~isfinite(lnw))=-inf;
    
  maxlnw=max(lnw);
    
  lnw=bsxfun(@minus,lnw,maxlnw);
    
  w=exp(lnw);
    
  w(not(isfinite(w)))=0;
    
  cumweights=cumsum(w);
    
  sumw=cumweights(end,:);
    
  logl=(log(sumw)+maxlnw)';
    
  if gpu_switch
        logl=gather(logl);
  end
    
  
  %%%%%%%%%%
  %resample%
  %%%%%%%%%%
  cumweights=bsxfun(@rdivide,cumweights,sumw);
  
  if gpu_switch
      randu=gpuArray.rand(R,Nparam);
  else
      randu=rand(R,Nparam);
  end
    
  v=(0:R-1)';
  
  if gpu_switch
      v=gpuArray(v);
  end
  
  u=(repmat(v,1,Nparam)+randu)/R;
  
  %row indices of u within cumweights 
  %values 1,... Nparticles
  %size R*Nparam
  bin=Resample_Mtx(u,cumweights);
  
  k=(0:Nparam-1);
  
  if gpu_switch
      k=gpuArray(k);
  end
  
  mtx_bin=repmat(k*Nparticles,R,1)+bin;
  