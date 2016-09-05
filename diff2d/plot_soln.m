function dummy = plot_soln(inp)

fid = fopen(inp);

m = fscanf(fid,'%d',1);
u = fscanf(fid,'%g',[m-1,m-1]);

fclose(fid);

% we don't bother with the Dirichlet bndry
x = (1/(m-2))*[0:m-2]';

surf(x,x,u);shading interp;
