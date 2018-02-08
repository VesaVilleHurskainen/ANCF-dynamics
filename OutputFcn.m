function status = OutputFcn( t, y, flag, data )
%OUTPUTFCN' Function that is called by the integrators in every timestep
%that generates an output.
if strcmp(flag, 'init')==1
else
    if ~isempty(t)
        tdisp=t(end);
		cpus = toc;
		cpum = floor(cpus/60);
		cpus = cpus - cpum*60;
		cpuh = floor(cpum/60);
		cpum = cpum - cpuh*60;
		clk = clock;
        display(['...t = ',num2str(tdisp),' s reached. Cumulative calculation time: ', num2str(cpuh),' h ', num2str(cpum), ' min ', num2str(cpus),' s. Printed at: ',num2str(clk(4)),':',num2str(clk(5)),':',num2str(floor(clk(6)))])

        % nx=data.nx;
        % lincon=data.lincon;
        % bc=data.bc;
        % ndof=data.ndof;
        % Case=data.Case;
        % elementID=data.elementID;
        % elemDim=data.elemDim;
        % ElemDofs=data.ElemDofs;
        % DofsAtNode=data.DofsAtNode;
        % nloc=data.nloc;
        % H=data.H;
        % W=data.W;
        % L=data.L;
        % resolution=data.resolution;
        % ee0=data.ee0;
        
        % figure(3)
        % eetemp=zeros(nx,1);
        % eetemp(lincon.dof)=ee0(lincon.dof)';
        % eetemp(bc)=y(1:ndof,1);
        % drawsys2d(Case,elementID,elemDim,ElemDofs,DofsAtNode,nloc,eetemp,H,W,L,0,resolution);
        % title(['t = ',num2str(t,'%5.2f')])
        % xlabel('x')
        % if Case==7
            % ylabel('r')
        % else
            % ylabel('y')
        % end
        % %axis equal
        % if ((Case==1) || (Case==2) || (Case == 3))
            % axis equal
            % xlim([-1.11,1.11])
            % ylim([-1.11,0.21])
        % elseif Case==7
            % xlim([-0.1,L+0.1])
            % ylim([-H/2-H,H/2+H])
        % end
        % grid on
        % drawnow expose update
        % drawnow;
    end
end
status = 0;
end

