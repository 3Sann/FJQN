function TH = AVMA(K,P,T)
    % the number of servers in service station
    rs = [K,P];
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
        P{i}{1}{1} = 1;
        Q{i}{1} = 0;
        L{i}{1} = 0;
        L_piao{i}{1} = 0;
    end

    % Preprocessing
    for i = 1:sta_inner
        temp_ES = (rs(i) - 1) / (rs(i) + 1) * ES(i) / rs(i) + 2 / (rs(i) + 1) / rs(i) * ES2(i) / (2 * ES(i));
        ES_rems(i) = temp_ES;
    end
    % Iteration 
    for t = 1:T
        for i = 1:sta_inner
            temp_ET = Q{i}{t}* ES_rems(i) + L_piao{i}{t} * ES(i) / rs(i) + ES(i);
            ET(i) = temp_ET;
        end
        TH = t / sum(vs .* ET);
        for i = 1:sta_inner
            b_up = min(rs(i) - 1, t); %求b的上线
            P{i}{t+1}{1} = 0;
            if b_up ~= 0                 
                 if b_up == 1
                   P{i}{t+1}{2} = ES(i) / 1 * vs(i) * TH * P{i}{1}{1};
                 else
                    for b = 1:b_up
                        temp_B = ES(i) / b * vs(i) * TH * P{i}{t}{b};
                        P{i}{t+1}{b+1} = temp_B;
                    end
                 end
            end
        end
        for i = 1:sta_inner
            if t < rs(i) 
                Q{i}{t+1} = 0;
            else
                if t == 1
                    Q{i}{t+1} = ES(i) / rs(i) * vs(i) * TH * (Q{i}{t}+ P{i}{t}{1});
                else
                    Q{i}{t+1} = ES(i) / rs(i) * vs(i) * TH * (Q{i}{t} + P{i}{t}{rs(i)});
                end
            end
            pro = sum(cell2mat(P{i}{t+1})) + Q{i}{t+1};
            if pro>1
                P{i}{t+1}{1} = 0;
            else
                P{i}{t+1}{1} = 1 - sum(cell2mat(P{i}{t+1})) - Q{i}{t+1};
            end
        end
        for i = 1:sta_inner
            if t < rs(i)
                L_piao{i}{t+1} = 0;
            else
                L_piao{i}{t+1} = ES(i) / rs(i) * vs(i) * TH * (Q{i}{t} + L_piao{i}{t});
            end
            L{i}{t+1} = TH * vs(i) * ET(i);
        end
    end
end
