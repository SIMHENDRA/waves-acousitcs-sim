%clear vars

clearvars testerjawn tester2jawn testerjawna tester2jawna hit count i jawn2 jawn3 testerjawnfreq jawn4 jawn5 tester2jawnfreq jawn5hitfreq jawn3hitfreq hitfreqpow1 hitfreqpow2 hitfreqAndpowstage1 hitfreqAndpowstage2 hitfreqAndpow diffs 




% find where peaks happen in fft of acous and s2
testerjawna = abs(fft(AcousSave,8192));
% testerjawn = zeros(size(testerjawna) - 1);
testerjawn = log10(testerjawna(2:end));

tester2jawna = abs(fft(S2,8192));
% tester2jawn = zeros(size(tester2jawna) - 1);
tester2jawn = log10(tester2jawna(2:end));


count = 1;
for i = 2:(length(testerjawn) - 1)
    
    if (testerjawn(i) > testerjawn(i-1) && testerjawn(i) > testerjawn(i+1))
        jawn2(count) = testerjawn(i);
        jawn3(count) = i;
        testerjawnfreq = [jawn2; jawn3];
        count = count + 1;
    end
end

count = 1;
for i = 2:(length(tester2jawn) - 1)
    
    if (tester2jawn(i) > tester2jawn(i-1) && tester2jawn(i) > tester2jawn(i+1))
        jawn4(count) = tester2jawn(i);
        jawn5(count) = i;
        tester2jawnfreq = [jawn4; jawn5];
        count = count + 1;
    end
end

% find freq where both acous and s2 fft's both

hit = 0;
for i = 1:length(jawn3)
    for b = 1:length(jawn5)
        if (abs(jawn5(b) - jawn3(i)) <= 3)
            hit = hit + 1;
            jawn5hitfreq(hit) = jawn5(b);
            jawn3hitfreq(hit) = jawn3(i);
            hitfreqpow1(hit) = testerjawnfreq(1,find(testerjawnfreq(2,:) == jawn3hitfreq(hit)));
            hitfreqpow2(hit) = tester2jawnfreq(1,find(tester2jawnfreq(2,:) == jawn5hitfreq(hit)));
            
%             hitfreqpow1(hit) = testerjawnfreq(1,find(testerjawnfreq(2,:) == jawn3hitfreq(hit)));
%             hitfreqpow2(hit) = tester2jawnfreq(1,find(tester2jawnfreq(2,:) == jawn5hitfreq(hit)));
            
        end
    end
end

hitfreqAndpowstage1 = [jawn3hitfreq; jawn5hitfreq];
hitfreqAndpowstage2 = [hitfreqAndpowstage1; hitfreqpow1];
hitfreqAndpow = [hitfreqAndpowstage2; hitfreqpow2];

diffs = abs(hitfreqAndpow(3,:) - hitfreqAndpow(4,:));

score = mean(diffs);

% freqHz2 = (0:1:length(tester2jawn)-1)*(8192/length(fft(S2)));
% freqHz = (0:1:length(testerjawn)-1)*(8192/length(fft(AcousSave)));
