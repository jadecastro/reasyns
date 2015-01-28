            if ~isempty(funnelIn{ii,trans{itrans}(1),k-1})
                count1 = count1+1;
                for j = 1:length(funnelIn{ii,trans{itrans}(1),k-1}.t)
                    tmp = inv(funnelIn{ii,trans{itrans}(1),k-1}.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    ellC11{count1}(j) = ellipsoid(funnelIn{ii,trans{itrans}(1),k-1}.x(j,:)',tmp*funnelIn{ii,trans{itrans}(1),k-1}.rho(j));
                end
                ellC1{count1}(j) = ellFunnel{ii,trans{itrans}(1),k-1};
                Vin1{count1} = ones(size(funnelIn{ii,trans{itrans}(1),k-1}.V)) - funnelIn{ii,trans{itrans}(1),k-1}.V;
            end
