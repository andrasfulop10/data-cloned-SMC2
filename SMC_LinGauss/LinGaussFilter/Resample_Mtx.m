function bin=Resample_Mtx(u,cumweights)
%find position of u in cumweights for each row, assumes u is ascending
%(stratified sampling)

[Nparticles,Nsample]=size(u);

jointmtx=[u; cumweights];

[PH,orig_position]=sort(jointmtx);

iscumweight=orig_position>Nparticles;

whichbin=cumsum(iscumweight)+1;

bin=reshape(whichbin(~iscumweight),Nparticles,Nsample);
