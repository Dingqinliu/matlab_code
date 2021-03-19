function nu=NPCR(P1,P2)
nu=zeros(1,1);
PR=P1(:,:,1);PR1=P2(:,:,1);
[M,N]=size(PR);
D=(PR~=PR1);
nu(1)=sum(sum(D))/(M*N)*100;
end