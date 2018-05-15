%Integrator function to numerically comput Satellite Dwarf Galaxy behavior
%in galalctic potential for use in ExtragalacitcBandStructuresv1.m

function [t,sdgn] = integrosdg(GM,xdist1,ydist1,zdist1,vxd1,vyd1,vzd1,Ndist,t1,t2,del)

[t,sdgn] = ode23(@sdgprim1,[t1,t2],[xdist1,ydist1,zdist1,vxd1,vyd1,vzd1]);

    function sdgprim = sdgprim1(ts,sdgn)
        sdgprim = zeros(6,1);
%         rd = zeros(1,Ndist+1);
%         Pd = zeros(1,Ndist+1);
%         for id = 2:Ndist+1
%             rd(id) = sqrt(xdist1(id)^2 + ydist1(id)^2 + zdist1(id)^2);
%             Pd(id) = sqrt(4*pi*pi*(rd(id)^3)/(GM(1)));
%         end %id loop
        
        %xs etc unnecessary for central potential, but allows us to move
        %location of potential if we wish.
        xs(1) = 0.0;
        ys(1) = 0.0;
        zs(1) = 0.0;
  
        sdgprim(1) = sdgn(4);
        sdgprim(2) = sdgn(5);
        sdgprim(3) = sdgn(6);

        sdgprim(4) = -GM(1)*(0.05+(sdgn(1)-xs(1))^2 + (sdgn(2)-ys(1))^2 + (sdgn(3)-zs(1))^2)^(-1-del)...
                    *(sdgn(1)-xs(1));
        sdgprim(5) = -GM(1)*(0.05+(sdgn(1)-xs(1))^2 + (sdgn(2)-ys(1))^2 + (sdgn(3)-zs(1))^2)^(-1-del)...
                    *(sdgn(2)-ys(1));
        sdgprim(6) = -GM(1)*(0.05+(sdgn(1)-xs(1))^2 + (sdgn(2)-ys(1))^2 + (sdgn(3)-zs(1))^2)^(-1-del)...
                    *(sdgn(3)-zs(1))  - 0.05*GM(1)*sign(sdgn(3));
                
    end %sdgprim
end %integrosdg

                
                
                
                