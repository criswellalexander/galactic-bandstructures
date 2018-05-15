function [t,uvn] = integro3d(del, Ndist, GM, t1, t2, x1, y1, vx1, vy1, z1, vz1,xdist2,ydist2,zdist2)

[t,uvn] = ode23(@uvprim1,[t1 t2],[x1 y1 vx1 vy1 z1 vz1]);

    function uvprim = uvprim1(tn,uvn)

        uvprim = zeros(6,1);    % a column vector
        ast(1) = 0.0;
        P(1) = 0.0;
        xs = xdist2;
        ys = ydist2;
        zs = zdist2;
        for id = 2:(Ndist+1) %Relative Distances
            dxs(id) = (uvn(1)-xs(id));
            dys(id) = (uvn(2)-ys(id));
            dzs(id) = (uvn(5)-zs(id));
            rxyz2(id) = dxs(id)^2 + dys(id)^2 +dzs(id)^2;
            rxyz(id) = sqrt(rxyz2(id));
        end %Rel dist loop
        
        uvprim(1) = uvn(3);
        uvprim(2) = uvn(4);
        uvprim(5) = uvn(6);
        uvprim(3) = -GM(1)*(0.01+(uvn(1)-xs(1))^2 + (uvn(2)-ys(1))^2 + (uvn(5)-zs(1))^2)^(-1-del)...
                *(uvn(1)-xs(1));
        uvprim(4) = -GM(1)*(0.05+(uvn(1)-xs(1))^2 + (uvn(2)-ys(1))^2 + (uvn(5)-zs(1))^2)^(-1-del)...
                *(uvn(2)-ys(1));
        uvprim(6) = -GM(1)*(0.05+(uvn(1)-xs(1))^2 + (uvn(2)-ys(1))^2 + (uvn(5)-zs(1))^2)^(-1-del)...
                *(uvn(5)-zs(1))  - 0.05*GM(1)*sign(uvn(5));

        for ig=2:(Ndist+1) %Disturber Gravity
            uvprim(3) = uvprim(3) - GM(ig)*((rxyz2(ig))^(-1.5))*dxs(ig);
             uvprim(4) = uvprim(4) - GM(ig)*((rxyz2(ig))^(-1.5))*dys(ig);
             uvprim(6) = uvprim(6) - GM(ig)*((rxyz2(ig))^(-1.5))*dzs(ig);
        end %Dist Grav loop
        clear xs ys zs dxs dys dzs rxyz2 rxyz
    end %uvprim
end

