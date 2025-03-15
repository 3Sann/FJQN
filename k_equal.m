function Kequal = k_equal(k)
    statusInFJ = nchoosek(4+k,4)*nchoosek(8+k,8) - nchoosek(4+k-1,4)*nchoosek(8+k-1,8);
    M = 12;
    statusAboutk = nchoosek(M - 1 + 1, M - 1);
    mindis = abs(statusInFJ - statusAboutk);
    Kequal = 1;
    for i = 2:1000
        statusAboutk = nchoosek(M - 1 + i, M - 1);
        dis = abs(statusInFJ - statusAboutk);
        if dis > mindis
            break
        end
        Kequal = i;
        if (Kequal - k)/k > 0.3
            break
        end
        mindis = dis;
    end
end
