function TH = OutMVA(K,P,T)
    vs = [1, 1, 1, 1];
    sta_inner = 4;
    for i = 1:sta_inner
        p{i}{1}{1} = 1;
    end
    for t=1:T
        for i =1:sta_inner
            ET{t}{i} = 0;
            for j = 0:t-1
                newT = k_equal(j+1,2);
                tempA = AVMA(K,P,newT);
                u = [tempA,tempA,(j+1)/6,(j+1)/6];
                ET{t}{i} = ET{t}{i} + (j+1)/u(i)*p{i}{t}{j+1};
            end
        end
        temp = 0;
        for i = 1:sta_inner
            temp = temp + vs(i)* ET{t}{i};
        end
        TH = t/temp;
        for i = 1:sta_inner
            temp = 0;
            for j=1:t
                newT = k_equal(j,2);
                tempA = AVMA(K,P,newT);
                u = [tempA,tempA,(j+1)/6,(j+1)/6];
                p{i}{t+1}{j+1} = vs(i)*TH/u(i)*p{i}{t}{j};
                temp = temp + p{i}{t+1}{j+1};
            end
            p{i}{t+1}{1} = 1 - temp;
        end
    end
end
