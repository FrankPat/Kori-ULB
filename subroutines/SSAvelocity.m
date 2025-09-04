function [uxssa,uyssa,beta2,eta,etaD,dudx,dudy,dvdx,dvdy,su,ubx,uby,ux,uy, ...
    damage,NumStabVel]= ...
    SSAvelocity(ctr,par,su,Hmx,Hmy,gradmx,gradmy,signx,signy,zeta, ...
    uxssa,uyssa,etaD,A3d,H,HB,B,stdB,Asf,A,MASK,glMASK,HAF,HAFmx,HAFmy,cnt, ...
    nodeu,nodev,MASKmx,MASKmy,bMASK,uxsia,uysia,udx,udy,ubx,uby,node,nodes, ...
    Mb,Melt,dtdx,dtdx2,VM,damage,shelftune)

% Kori-ULB
% Iterative solution to the SSA velocity (with SSA, hybrid and DIVA model)
% ctr.SSA=1: uxssa, uyssa = SSA velocity
% ctr.SSA=2: uxssa, uyssa = basal SSA velocity (without deformational)
% ctr.SSA=3: uxssa, uyssa = total vertically integrated velocity

    taudx=par.rho*par.g*Hmx.*sqrt(gradmx).*signx;
    taudy=par.rho*par.g*Hmy.*sqrt(gradmy).*signy;
    if ctr.uSSAexist==1 || cnt>1
        ussa=vec2h(uxssa,uyssa);    %VL: ussa on h-grid
    else
        ussa=zeros(ctr.imax,ctr.jmax)+0.1;
    end
    ussa=max(ussa,1e-3);
    if par.ShelfPinning==1 && ctr.stdBexist==1 && ctr.inverse==0 
        fg=max(0,1-(HB-B)./stdB); % subgrid pinning points in ice shelf
    else
        fg=1;
    end
    if ctr.u0>1e10
        beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
    else
        beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).* ...
            Asf/ctr.u0).^(-1/ctr.m);
    end
    if ctr.SSA==3
        % Integration factor.
        [F2,~,~]=Fint(ctr,etaD,H,zeta);
        % Effective beta from F2 correction in DIVA model.
        beta2=beta2./(1.0+beta2.*F2);
        beta2(beta2>1e8)=1./F2(beta2>1e8); % frozen conditions
    end
    beta2=min(beta2,1e8);
    beta2(MASK==0)=0;
    betax=0.5*(beta2+circshift(beta2,[0 -1]));
    betay=0.5*(beta2+circshift(beta2,[-1 0]));

    if ctr.mismip>=1
        betax(:,1)=betax(:,2); % symmetric divide
        betax(1,:)=betax(3,:); % symmetry axis
        betax(ctr.imax,:)=betax(ctr.imax-2,:); % periodic BC
        betax(:,ctr.jmax)=0; % ocean
        betay(:,1)=betay(:,3); % symmetric divide
        betay(1,:)=betay(2,:); % symmetry axis
        betay(ctr.imax,:)=betay(ctr.imax-1,:); % periodic BC
        if ctr.mismip==2 % Thule setup
            betax(ctr.imax,:)=0;
            betay(ctr.imax,:)=0;
        end
    end

    % Viscosity on edge (or when viscosity cannot be calculated)
    if ctr.SSA==2 % hydrid model
        udx(HAFmx<0)=0; % no deformation for floating cells
        udy(HAFmy<0)=0;
    else % pure SSA
        udx=zeros(ctr.imax,ctr.jmax);
        udy=zeros(ctr.imax,ctr.jmax);
    end
    for ll=1:par.visciter % iteration over effective viscosity
        if ctr.SSA<3
            % Effective viscosity for SSA and hybrid model
            [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,MASK, ...
                glMASK,shelftune,damage,ctr);
        else
            % Effective viscosity for DIVA solver
            [eta,etaD,dudx,dvdy,dudy,dvdx]=EffViscDIVA(A3d,betax,betay,ubx,uby, ...
                etaD,H,damage,uxssa,uyssa,zeta,MASK,glMASK,shelftune,ctr,par);
        end
        if ctr.damage==1 && (ctr.damexist==1 || cnt>1)
            if ll==1
                dtr=zeros(ctr.imax,ctr.jmax);
                ThinComp=zeros(ctr.imax,ctr.jmax);
                if ctr.TRdam==1
                    % compute advected damage field from previous timestep (dtr)
                    if ctr.THdam==1 
                        % Thinning component (very sensitive)
                        ThinComp=ThinningComponent(ctr,par,dudx,dvdy, ...
                            dudy,dvdx,eta,H,damage,bMASK,glMASK);
                    end
                    SuM=Mb; BoM=Melt;
                    if ctr.SFdam==0
                        SuM=SuM*0;
                    end
                    if ctr.BSdam==0
                        BoM=BoM*0;
                    end
                    dtr=TransportDamage(node,nodes,damage,SuM, ...
                        BoM,ThinComp,H,glMASK,dtdx,dtdx2, ...
                        uxssa,uyssa,ctr,cnt,bMASK,VM,par);
                end
                % compute surface damage
                ds=SurfaceDamage(ctr,par,dudx,dvdy,dudy,dvdx,eta,H);
                % compute basal damage (and Kachuck term, necessary for transport)
                db=BasalDamage(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,HAF);
                dloc=max(0,min(db+ds,H.*par.damlim)); % local damage
                % final damage field taken as the maximum of dloc and dtr
                % allows for adverction of damage into regions that
                % would not initiate damage total damage is limited to damlim
                damage=min(H.*par.damlim,max(dloc,dtr));
                % ensure no damage in the open ocean
                damage(glMASK==6)=0;
                if ctr.basin==1
                    % ensure no outise of basin boundary
                    damage(bMASK==1)=0;
                end
                % scaling of viscosity
                scale_eta=1;
            end
        else
            scale_eta=1;
            if ctr.damexist==0
                % If damage exists, this will keep damage constant, else
                % initializes it to zero damage
                damage=zeros(ctr.imax,ctr.jmax);
            end
        end
        eta=eta.*scale_eta;
        
        [uxs1,uys1,su,flagU,relresU,iterU]=SparseSolverSSA(nodeu,nodev, ...
            su,MASKmx,MASKmy,bMASK,glMASK,H,eta,betax,betay,uxssa,uyssa, ...
            uxsia,uysia,udx,udy,taudx,taudy,ctr,par);
        duxs=sqrt((uxs1-uxssa).^2+(uys1-uyssa).^2);
        duxs(isnan(duxs))=0;
        uxssa=uxs1;
        uyssa=uys1;
        limit=sum(duxs(:))/(ctr.imax*ctr.jmax);
        if ctr.SSA==3 % DIVA solver
            % Integration factor.
            [F2,F2x,F2y]=Fint(ctr,etaD,H,zeta);
            % Eq. 32, Lipscomb et al. (2019).
            ubx=uxssa./(1.0+betax.*F2x);
            uby=uyssa./(1.0+betay.*F2y);
        end        
        %---------iterative beta---------
        if cnt<=ctr.BetaIter
            ussa=vec2h(uxssa,uyssa); %VL: ussa on h-grid
            if ctr.SSA==3
                ussa=vec2h(ubx,uby); % DIVA: initialize with basal velocity
            end
            if ctr.u0>1e10
                beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
            else
                beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).*Asf ...
                    /ctr.u0).^(-1/ctr.m);
            end
            if ctr.SSA==3
                % Effective beta from F2 correction in DIVA model.
                beta2=beta2./(1.0+beta2.*F2);
                beta2(beta2>1e8)=1./F2(beta2>1e8); % frozen conditions
            end
            beta2=min(beta2,1e8);
            beta2(MASK==0)=0;
            betax=0.5*(beta2+circshift(beta2,[0 -1]));
            betay=0.5*(beta2+circshift(beta2,[-1 0]));
            if ctr.mismip>=1
                betax(:,1)=betax(:,2); % symmetric divide
                betax(1,:)=betax(3,:); % symmetry axis
                betax(ctr.imax,:)=betax(ctr.imax-2,:); % periodic BC
                betax(:,ctr.jmax)=0; % ocean
                betay(:,1)=betay(:,3); % symmetric divide
                betay(1,:)=betay(2,:); % symmetry axis
                betay(ctr.imax,:)=betay(ctr.imax-1,:); % periodic BC
                if ctr.mismip==2 % Thule setup
                    betax(ctr.imax,:)=0;
                    betay(ctr.imax,:)=0;
                end
            end
        end
        %--------------------------------
        if limit<par.visctol % Limit on convergence
            break;
        end
    end
    if ctr.NumCheck==1
        NumStabVel=[ll,limit,relresU,iterU,flagU];
    else
        NumStabVel=false;
    end
    uxssa(Hmx==0 & HAFmx>0)=0; % only for grounded ice sheet with H=0
    uyssa(Hmy==0 & HAFmy>0)=0;
    uxssa=min(max(-par.maxspeed,uxssa),par.maxspeed);
    uyssa=min(max(-par.maxspeed,uyssa),par.maxspeed);
    ux=uxssa;
    uy=uyssa;
    if ctr.SSAdiffus==2 % SSA as basal velocity only
        ux=ux+udx;
        uy=uy+udy;
    end
    if ctr.SSA==1 % basal sliding velocity SSA
        ubx=ux;
        uby=uy;
    end
    if ctr.SSA==2
        ubx=ux-udx;
        uby=uy-udy;
    end
    
    % return beta2 value instead of effective beta2 for DIVA model
    if ctr.SSA==3
        beta2=beta2./(1-beta2.*F2);
        beta2(beta2>1e8)=1./F2(beta2>1e8);
        beta2=min(beta2,1e8);
        beta2(MASK==0)=0;
    end
    
end


