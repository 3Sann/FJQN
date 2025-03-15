function TH = AVMA(K)
    % the number of servers in service station
    rs = [K,6];
    % the arriveal rate of service station
    vs = [1, 1];
    % the first moment, second moment of the service station
    ES = [13,311,11,1247];
    ES2 = [8/397,120/986];
    % the number of stations
    sta_inner = 2;
    TH = 0;
    
    % Initialize
    for i = 1:sta_inner
        temp_p0 = 1;
        temp_q0 = 0;
        temp_l0 = 0;
        temp_l0_piao = 0;
        P0{i} = temp_p0;
        Q0{i} = temp_q0;
        L0{i} = temp_l0;
        L0_piao{i} = temp_l0_piao;
    end
    Q{1} = Q0;
    L{1} = L0;
    P{1} = P0;
    L_piao{1} = L0_piao;

    % Preprocessing
    for i = 1:sta_inner
        temp_ES = (rs(i) - 1) / (rs(i) + 1) * ES(i) / rs(i) + 2 / (rs(i) + 1) / rs(i) * ES2(i) / (2 * ES(i));
        ES_rems{i} = temp_ES;
    end
    index = 2;
    % Iteration 
    for k = 1:K
        for i = 1:sta_inner
            temp_ET = Q{index-1}{i} * ES_rems{i} + L_piao{index-1}{i} * ES(i) / rs(i) + ES(i);
            ET{i} = temp_ET;
        end
        TH = k / sum(vs .* cell2mat(ET));
        for i = 1:sta_inner
            b_up = min(rs(i) - 1, k); %求b的上线
            P{index}{i}{1} = 0;
            if b_up ~= 0                 
                 if b_up == 1
                   P{index}{i}{2} = ES(i) / 1 * TH * P{1}{i};
                 else
                    for b = 1:b_up
                        temp_B = ES(i) / b * TH * P{index-1}{i}{b};
                        P{index}{i}{b+1} = temp_B;
                    end
                 end
            end
        end
        for i = 1:sta_inner
            if k < rs(i)
                Q{index}{i} = 0;
            else
                if k == 1
                     Q{index}{i} = ES(i) / rs(i) * vs(i) * TH * (Q{index-1}{i} + P{index-1}{i});
                else
                    Q{index}{i} = ES(i) / rs(i) * vs(i) * TH * (Q{index-1}{i} + P{index-1}{i}{rs(i)});
                end
            end
            P{index}{i}{1} = 1 - sum(cell2mat(P{index}{i})) - Q{index}{i};
        end
        for i = 1:sta_inner
            if k < rs(i)
                L_piao{index}{i} = 0;
            else
                L_piao{index}{i} = ES(i) / rs(i) * vs(i) * TH * (Q{index-1}{i} +L_piao{index-1}{i});
            end
            L{index}{i} = TH * vs(i) * ET{i};
        end
        index= index +1;
    end
end
