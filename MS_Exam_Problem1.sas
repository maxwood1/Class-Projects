*Applied STAT MS Exam ANOVA Model Designs for Problem 1;
*Max Woodbury;


*create fake data for each design;

*single measurement WITHOUT patient variability;
data singlenovar;
	call streaminit(200);
	do isim = 1 to 1000; *number of simulations;
	trt=1;
	do patient = 1 to 40;
		if trt > 12 then trt=1;
		if catheter > 1 then catheter=0;
		do organ = 1 to 2;
			do sensor = 1 to 3;
				y=rand('normal', 50, 10); *random noise;
				if catheter = 1 then y+4; 
				if organ = 1 then y+4;
				if sensor = 1 then y+5;
				if sensor = 2 then y+2.5;
				output;
				trt+1;
			end;
		end;
		catheter+1;
	end;
	end;
run;

*single measurement WITH patient variability;
data singlevar;
	call streaminit(200);
	do isim = 1 to 1000; *number of simulations;
	trt=1;
	do patient = 1 to 40;
		u = floor(rand('uniform')*3);
		if trt > 12 then trt=1;
		if catheter > 1 then catheter=0;
		do organ = 1 to 2;
			do sensor = 1 to 3;
				y=rand('normal', 50, 10); *random noise;
				if u = 0 then y=y-10; *patient variability;
				if u = 1 then y=y+10;
				if catheter = 1 then y+4;
				if organ = 1 then y=y+4;
				if sensor = 1 then y+5;
				if sensor = 2 then y+2.5;
				output;
				trt+1;
			end;
		end;
		catheter+1;
	end;
	end;
run;

*multiple measurements per sensor WITH patient variability;
data multi;
	call streaminit(200);
	do isim = 1 to 1000; *number of simulations;
	trt=1;
	do patient = 1 to 40;
		u = floor(rand('uniform')*3);
		if trt > 12 then trt=1;
		if catheter > 1 then catheter=0;
		do organ = 1 to 2;
			do sensor = 1 to 3;
				do rep = 1 to 2;
					y=rand('normal', 50, 10); *random noise;
					if u = 0 then y=y-10; *patient variability;
					if u = 1 then y=y+10;
					if catheter = 1 then y+4;
					if organ = 1 then y=y+4;
					if sensor = 1 then y+5;
					if sensor = 2 then y+2.5;
					output;
				end;
				trt+1;
			end;
		end;
		catheter+1;
	end;
	end;
run;



*Do power simulations for each design;

*Turn SAS output off for simulations;
ods graphics off;
ods exclude all;
ods results off;
options nonotes;

proc glm data=singlenovar;
	class catheter organ sensor;
	model y = catheter|organ|sensor / e3;
	estimate 'C1 vs C2' catheter 1 -1;
	contrast 'C1 vs C2' catheter 1 -1;
	estimate 'O1 vs O2' organ 1 -1;
	estimate 'S1 VS S2' sensor 1 -1 0;
    estimate 'S1 VS S3' sensor 1 0 -1;
    estimate 'S2 VS S3' sensor 0 1 -1;
	by isim;
	ods output ModelANOVA=anova;
run;
quit;

data power;
	set anova;
	where source = "catheter" or source = "organ" or source = "sensor";
	alpha = 0.05;
	fcrit = finv(1- alpha, df, 228, 0);
	noncent = df*fvalue;
	power = 1 - probf(fcrit, df, 228, noncent);
	if noncent > 3000 then power = 1.0;
run;

proc sort data=power;
	by source;
run;

*Turn SAS output back on;
ods graphics on;
ods exclude none;
ods results on;
options notes;

proc means data=power;
 var power;
 by source;
run;


*Turn SAS output off for simulations;
ods graphics off;
ods exclude all;
ods results off;
options nonotes;

proc glm data=singlevar;
	class patient catheter organ sensor;
	model y = catheter|organ|sensor patient(catheter) / e3;
	random patient(catheter);
	estimate 'C1 vs C2' catheter 1 -1;
	estimate 'O1 vs O2' organ 1 -1;
	estimate 'S1 VS S2' sensor 1 -1 0;
    estimate 'S1 VS S3' sensor 1 0 -1;
    estimate 'S2 VS S3' sensor 0 1 -1;
	by isim;
	ods output ModelANOVA=anova;
run;
quit;

data power1;
	set anova;
	where source = "catheter" or source = "organ" or source = "sensor";
	alpha = 0.05;
	fcrit = finv(1- alpha, df, 190, 0);
	noncent = df*fvalue;
	power = 1 - probf(fcrit, df, 190, noncent);
	if noncent > 3000 then power = 1.0;
run;

proc sort data=power1;
	by source;
run;

*Turn SAS output back on;
ods graphics on;
ods exclude none;
ods results on;
options notes;

proc means data=power1;
 var power;
 by source;
run;



*Turn SAS output off for simulations;
ods graphics off;
ods exclude all;
ods results off;
options nonotes;

proc glm data=multi;
	class rep catheter organ sensor patient;
	model y = catheter|organ|sensor rep patient(catheter) / e3;
	random rep patient(catheter);
	estimate 'C1 vs C2' catheter 1 -1;
	estimate 'O1 vs O2' organ 1 -1;
	estimate 'S1 VS S2' sensor 1 -1 0;
    estimate 'S1 VS S3' sensor 1 0 -1;
    estimate 'S2 VS S3' sensor 0 1 -1;
	by isim;
	ods output ModelANOVA=anova;
run;
quit;

data power2;
	set anova;
	where source = "catheter" or source = "organ" or source = "sensor";
	alpha = 0.05;
	fcrit = finv(1- alpha, df, 429, 0);
	noncent = df*fvalue;
	power = 1 - probf(fcrit, df, 429, noncent);
	if noncent > 3000 then power = 1.0;
run;

proc sort data=power2;
	by source;
run;

*Turn SAS output back on;
ods graphics on;
ods exclude none;
ods results on;
options notes;

proc means data=power2;
 var power;
 by source;
run;
