function bnd = decouplesurf(bnd)
for ii = 1:length(bnd)-1
  % Despite what the instructions for surfboolean says, surfaces should be ordered from inside-out!!
  [newnode, newelem] = surfboolean(bnd(ii+1).pos,bnd(ii+1).tri,'decouple',bnd(ii).pos,bnd(ii).tri);
  bnd(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(bnd(ii+1).pos,1);
  bnd(ii+1).pos = newnode(newnode(:,4)==2,1:3);
end % for
end % function