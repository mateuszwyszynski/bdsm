
@----------------------------------------------------------------------@
@                                                                      @
@                     OPTIMIZATION PROCEDURES                          @
@                                                                      @
@----------------------------------------------------------------------@


proc (4) = optmum(fnct,x0);
    local x,f0,g,retcode;
    local Lopfhess,Lopitdta,LLoutput;
    local _opalgr, _opdelta,_opdfct,_opditer,_opgdmd,_opgdprc,_opgrdh,_opgtol,_ophsprc,
      _opkey,_opmbkst,_opmdmth,_opmiter,_opmtime,_opmxtry,_opparnm,_oprteps,_opshess,
      _opstep,_opstmth,_opusrch,_opusrgd,_opusrhs,__title,_opitdta,_opfhess;


    _opalgr = 2;    /* optimization algorithm */
    _opparnm = 0;           /* parameter names */
    _opstep = 2;            /* selects type of step length */
    _opshess = 0;           /* selects starting hessian */
    _opmbkst = 10;          /* # of backsteps in computing steplength  */
    _opgtol = 1e-5;         /* convergence tolerance for gradient */
    _ophsprc = 0;           /* procedure to compute hessian */
    _opgdprc = 0;           /* procedure to compute gradient */
    _opgdmd = 0;            /* numerical gradient method */
    _opditer = 20;          /* # iters to switch algorithms for _opmdmth  */
    _opdfct = 0.001;        /* % change in function for _opmdmth */
    _opmiter = 1e+5;        /* maximum number of iterations */
    _opmtime = 1e+5;
    _opitdta = { 0,0,0 };
    _oprteps = .01;
    _opusrch = 0;
    _opdelta = .1;
    _opmxtry = 100;
    _opfhess = 0;
    _opusrgd = 0;
    _opusrhs = 0;
    _opstmth = "";
    _opmdmth = "";
    _opkey = 1;
    _opgrdh = 0;
    __title = "";

#IFUNIX
    LLoutput = __output /= 0;
#else
    LLoutput = __output;
#endif

    { x,f0,g,retcode,Lopfhess,Lopitdta } = _optmum(fnct,x0,
      _opalgr,
      _opdelta,
      _opdfct,
      _opditer,
      _opgdmd,
      _opgdprc,
      _opgrdh,
      _opgtol,
      _ophsprc,
      _opkey,
      _opmbkst,
      _opmdmth,
      _opmiter,
      _opmtime,
      _opmxtry,
      _opparnm,
      _oprteps,
      _opshess,
      _opstep,
      _opstmth,
      _opusrch,
      _opusrgd,
      _opusrhs,
      LLoutput,
      __title );

    _opfhess = Lopfhess;
    _opitdta = Lopitdta;

    retp(x,f0,g,retcode);
endp;


#include gauss.ext


proc (6) = _optmum(fnct,x0,
    Lopalgr,
    Lopdelta,
    Lopdfct,
    Lopditer,
    Lopgdmd,
    Lopgdprc,
    Lopgrdh,
    Lopgtol,
    Lophsprc,
    Lopkey,
    Lopmbkst,
    Lopmdmth,
    Lopmiter,
    Lopmtime,
    Lopmxtry,
    Lopparnm,
    Loprteps,
    Lopshess,
    Lopstep,
    Lopstmth,
    Lopusrch,
    Lopusrgd,
    Lopusrhs,
    LLoutput,
    LLtitle
    );

    /* ------- LOCALS ----------- */
    local x,g,s,h,iter,ky,old,vof,d,dfct,f0,parnms,bksteps,dx,
        smallval,relgrad,algrm,stepm,pg,k0,k1,k2,lr,lf,ll,np,rteps,
        tstart,oldt,w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,chiter,h0,
        ttime,isctu,mask,fmt,omat,ret,hsmeth,beta,vof1,vofs,ds,it,oldfmt,
        tmpop,e1,e2;
    local midalg,midstp,midhs,hessm,c0,c1,skip;
    local Lopfhess,Lopitdta;
    local fnct:proc;
    clear bksteps,s,dfct,f0,h,chiter,iter,skip,Lopfhess,beta;
    x0 = vec(x0);
    if Lopusrgd /= 0;
       Lopgdmd = 2;
    endif;
    dx = 1;
    isctu = 1;
    pg = 1;
    smallval = 1e-15;
    if LLoutput == 2;
        call csrtype(0);
    endif;
    let w1  = 2  1  2 80 0   7;
    let w2  = 5  1 24 80 0   7;
    let w3  = 1  1  1 80 0 112;
    let w4  = 3  1  4 80 0 112;
    let w5  = 6  1 25 80 0   7;
    let w9  = 2  1  2  9 0   7;
    let w10 = 6 37  6 80 0   7;
    w0 = 196*ones(1,24);
    w6 = chrs(218~w0~194~w0~194~w0~191);
    w7 = chrs(179~(32*ones(1,24))~179~(32*ones(1,24))~179~(32*ones(1,24))~179);
    w8 = chrs(192~w0~193~w0~193~w0~217);
    w0 = "" $+ chrs(32*ones(40-floor(strlen(LLtitle)/2),1)) $+ LLtitle;

    algrm = "STEEP"|"BFGS"|"BFGS-SC"|"DFP"|"NEWTON"|"PRCG"|"NR";
    stepm = "1.0"|"STEPBT"|"HALF"|"BRENT"|"1"|"ONE"|"GOLDEN";
    hessm = "NOHESS"|"HESS";
    old = ndpcntrl(0,0);
    call ndpcntrl(1,1);
    if LLoutput == 2;
        cls;
    endif;
    clear midalg,midstp,midhs;
    if Lopstmth $/= "";
        Lopalgr = _strsch(Lopstmth,algrm);
        Lopstep = _strsch(Lopstmth,stepm);
        Lopshess = _strsch(Lopstmth,hessm);
    endif;

    gosub chk1(Lopalgr,2);
        pop Lopalgr;
    gosub chk2(Lopstep,2);
        pop Lopstep;
    gosub chk4(Lopshess);
        pop Lopshess;

    if Lopmdmth $/= "";
        gosub chk1(_strsch(Lopmdmth,algrm),Lopalgr);
            pop midalg;
        gosub chk2(_strsch(Lopmdmth,stepm),Lopstep);
            pop midstp;
        gosub chk3(_strsch(Lopmdmth,hessm),Lopshess);
            pop midhs;
    endif;
    if Lopparnm $/= 0 and rows(Lopparnm) /= rows(x0);
        if not trapchk(4);
            if LLoutput == 2;
                locate 2,1;
            endif;
            errorlog "vector of parameter labels does not conform to"\
                     "vector of starting values";
        endif;
        parnms = 0;
    else;
        parnms = Lopparnm;
    endif;
    if parnms $== 0;
        let mask[1,3] = 1 1 1;
        let fmt[3,3] = "lf " 8 0 "lf" 18 4 "lf" 18 4;
    else;
        let mask[1,3] = 0 1 1;
        let fmt[3,3] = "s " 8 8 "lf" 18 4 "lf" 18 4;
    endif;

/****************************************************************************/
/*                     BEGIN OPTIMIZATION                                   */
/****************************************************************************/
    tstart = date;
    np = rows(x0);          /* Number of parameters to estimate */

    if LLoutput == 2;
        gosub BAR;
    endif;
    ttime = date;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    vof = fnct(x0);       /* Initial function value */
    if scalmiss(vof) or
       vof == __INFp or
       vof == __INFn or
       vof == __INDEFn or
       vof == __INDEFp;
        if not trapchk(4);
            if LLoutput == 2;
                locate 2,1;
            endif;
            errorlog "ERROR:  function cannot be computed at initial"\
                     " parameter values";
        endif;
        ret = error(3);
        x = x0;
        g = miss(zeros(rows(x),1),0);
        goto A98;
    endif;
    if LLoutput == 2;
        scroll w9;
    endif;

    if LLoutput == 2;
        locate 2,15;
        printdos " Gradient";
    endif;
    g = _deriv2(x0,1,&fnct,Lopgdmd,Lopgdprc,Lophsprc,LLoutput,
                Lopusrgd,Lopusrhs,Lopgrdh)';
    if scalmiss(g);
        if not trapchk(4);
            errorlog "gradient function failed at initial values";
        endif;
        ret = error(4);
        x = x0;
        goto A98;
    endif;
    if not(rows(x0) == 1 and rows(x0) == 1);
        if rows(g) == 1 and cols(g) /= 1;
            if not trapchk(4);
                errorlog "The gradient function has returned a column vector"\
                         " rather than the required row vector";
            endif;
            ret = error(9);
            x = x0;
            goto A98;
        endif;
    endif;
    if (rows(g)/=rows(x0));
        if not trapchk(4);
            errorlog "The number of elements in the gradient function";
            errorlog "is inconsistent with the number of starting values";
        endif;
        ret = error(8);
        x = x0;
        goto A98;
    endif;
    relgrad = (abs(g).*maxc(abs(x0)'|ones(1,rows(x0))))/maxc(abs(vof)|1);
    if abs(g) < smallval or relgrad < Lopgtol or Lopmiter == 0;
        x = x0;
        if Lopmiter == 0;
            ret = error(2);
        else;
            ret = error(0);
        endif;
        goto A98;
    endif;

    if rows(Lopshess) == 1;
        if LLoutput == 2;
            scroll w1;
        endif;
        if Lopshess == 0 and Lopalgr >= 2 and Lopalgr <= 4;
            if LLoutput == 2;
                locate 2,1;
                printdos " H set to identity matrix";
            endif;
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        elseif Lopalgr >= 2 and Lopalgr <= 5;
            gosub hssn;
               pop h;
        elseif Lopalgr == 6;
            h = 1;
        endif;
    elseif Lopalgr <= 4;
        h = chol(Lopshess);
    elseif Lopalgr == 5;
        h = Lopshess;
    endif;

    if np gt 48 and LLoutput == 2;
        locate 25,15;
        printdos "\27[7m";
        print "  <PgDn>, <PgUp>  page parameters and gradient ";
        printdos "\27[0m";
    endif;
    ttime = date;

A0:

/* ********* Start of iteration loop ********** */
    iter = iter + 1;
    f0 = vof;
    if LLoutput == 2;
        scroll w1;
        gosub BAR;
        gosub PARBOX;
        printdos "\27[7m";
        locate 3,29;
        printdos ftos(vof,"%-*.*lf",10,5);
        locate 4,13;
        printdos ftos(ethsec(tstart,date)/100,"%-*.*lf",4,2);
        locate 3,52;
        printdos algrm[Lopalgr];
        locate 3,66;
        printdos stepm[Lopstep];
        locate 4,70;
        printdos ftos(s,"%-*.*lf",5,3);
        if iter > 1;
            locate 4,27;
            printdos ftos(dfct,"%*.*lE",10,2);
        endif;
        locate 4,52;
        printdos ftos(bksteps,"%-*.*lf",1,0);
        locate 3,7;
        printdos ftos(iter,"%-*.*lf",1,0);
        printdos "\27[0m";

        k0 = lf;
        k2 = 0;
        do until k2 == lr;
            k2 = k2+1;
            k1 = 0;
            do until k1 == 3;
                k1 = k1+1;
                locate 7+k2,(k1-1)*25+2;
                if parnms == 0;
                    printdos ftos(k0,"%*.*lf",3,0);
                    printdos ftos(x0[k0],"%*.*lf",9,4);
                    printdos ftos(relgrad[k0],"%*.*lf",9,5);
                else;
                    printdos parnms[k0];
                    printdos " ";
                    printdos ftos(x0[k0],"%*.*lf",6,3);
                    printdos ftos(relgrad[k0],"%*.*lf",7,4);
                endif;
                k0 = k0+1;
                if k0 gt np;
                    goto A1;
                endif;
            endo;
        endo;
    elseif LLoutput == 1;
        print;
        print "============================================================"\
            "====================";
        oldfmt = sysstate(19,0);
        format /ld 4,0;
        print "   iteration: " iter;
        format 6,6;
        print "   algorithm: " $algrm[Lopalgr];;
        print "      step method: " $stepm[Lopstep];
        format /ld 10,5;
        print "   function: " vof;;
        print "   step length: " s;;
        format 3,0;
        print "   backsteps: " bksteps;
        print "------------------------------------------------------------"\
            "--------------------";
        print "   param.      param. value     relative grad.";
        call sysstate(19,oldfmt);

        if parnms $== 0;
            omat = seqa(1,1,np)~x0~relgrad;
        else;
            omat = parnms~x0~relgrad;
        endif;
        call printfm(omat,mask,fmt);
    endif;
A1:

    tstart = date;

/* =================== DIAGNOSTIC =================== */

/*
**   if you wish to inspect one of the following:
**
**    coefficient vector             x0
**    value of the function          vof
**    gradient                       g
**    Hessian                        h
**
**   then uncomment the appropriate line below:
*/

/*    print "coefficients " x0;   */
/*    print "function " vof;      */
/*    print "gradient " g;        */

/*
      if Lopalgr > 1 and Lopalgr < 5;
         print "Hessian estimate " h'h;
      elseif Lopalgr == 5;
         print "Hessian " h;
      endif;
*/


/*
**   if your program is crashing and you wish to retrieve one of
**   these arrays at the point of the crash, uncomment the appropriate
**   line below:
*/

/*     clearg x0_sav;  x0_sav = x0;     */
/*     clearg vof_sav; vof_sav = vof;   */
/*     clearg g_sav;   g_sav = g;       */

/*
      clearg h_sav;
      if Lopalgr > 1 and Lopalgr < 5;
         h_sav = h'h;
      elseif Lopalgr == 5;
         h_sav = h;
      endif;
*/

/* =================== end DIAGNOSTIC =================== */

    if Lopalgr == 1;
        d = -g;
    elseif Lopalgr >= 2 and Lopalgr <= 4;
        d = -cholsol(g,h);
    elseif Lopalgr == 5;
        oldt = trapchk(1);
        trap 1,1;
        d = -solpd(g,h);
        trap oldt,1;
        if scalmiss(d);
            if not trapchk(4);
                if LLoutput == 2;
                    locate 2,15;
                endif;
                errorlog "Newton iteration failed";
            endif;
           if Lopdelta /= 0;
                { e1,e2 } = eigrs2(h);
                e1 = 1/(e1 + Lopdelta - minc(e1));
                d = -(e2*diagrv(eye(rows(e1)),e1)*e2')*g;
            else;
                d = -g;
            endif;
        endif;
    elseif Lopalgr == 6;
        d = -g + beta;
    endif;

    { s,bksteps } = _stepl2(g,vof,x0,d,&fnct,Lopstep,Lopusrch,Lopmbkst,
                                                        LLoutput,Lopmxtry);
    if not scalmiss(s);
        if LLoutput == 2;
            locate 2,1;
            printdos " Function";
        endif;
        vof1 = fnct(x0 + s*d);
        if vof1 == __INFp or
           vof1 == __INFn or
           vof1 == __INDEFn or
           vof1 == __INDEFp;
           s = error(0);
        endif;
    else;
       vof1 = 0;
    endif;
    if scalmiss(s) or vof1 > vof;
        if Loprteps;
            s = 1;
            vof1 = vof + 1;
            vofs = 1e200;
            ds = 1;
            if LLoutput == 2;
                locate 2,15;
                printdos " Random Search";
            endif;
            it = 1;
            do while vof1 >= vof and it < Lopmxtry;
                if LLoutput == 2;
                    locate 2,30;
                    printdos ftos(it,"%-*.*lf",1,0);
                endif;
                rteps = 10^trunc(log(meanc(abs(g)))) * Loprteps;
                d = rteps*(2*rndu(rows(d),1)-1).*x0;
                vof1 = fnct(x0+d);
                if scalmiss(vof1) or
                    vof1 == __INFp or
                    vof1 == __INFn or
                    vof1 == __INDEFn or
                    vof1 == __INDEFp;
                     vof1 = vofs;
                     d = ds;
                elseif vof1 < vofs;
                   vofs = vof1;
                   ds = d;
                endif;
                it = it + 1;
            endo;
            if it >= Lopmxtry;
                vof1 = vofs;
                d = ds;
            endif;
            if vofs == 1e200;
               ret = error(3);
               x = x0;
               goto A98;
            endif;
            if LLoutput == 2;
                 scroll w1;
            endif;
            if Lopalgr > 1 and Lopalgr < 5;
                h = eye(np)*maxc(sqrt(abs(vof1))|1);
                isctu = 0;
            elseif Lopalgr == 5;
                isctu = 1;
            elseif Lopalgr == 6;
                beta = 0;
                isctu = 0;
            endif;
        else;
            if not trapchk(4);
                if scalerr(s) == 6;
                    errorlog "step length calculation failed";
                elseif scalerr(s) == 3;
                    errorlog "function calculation failed";
                endif;
            endif;
            if LLoutput == 2;
                cls;
            endif;
            x = x0;
            ret = s;
            goto A98;
        endif;
    endif;
    dx = s*d;
    x = x0 + dx;
    x0 = x;
    vof = vof1;
    if LLoutput == 2;
        scroll w9;
    endif;
    if isctu;
        { g,h,beta } = _sctu2(x,vof,smallval,g,h,dx,d,s,&fnct,Lopalgr,
                 LLoutput,Lopgdmd,Lopgdprc,Lophsprc,Lopusrgd,Lopusrhs,Lopgrdh);
        if scalmiss(g);
            ret = error(4);
            goto A98;
        elseif scalmiss(h);
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        endif;
    else;
        g = _deriv2(x,1,&fnct,Lopgdmd,Lopgdprc,Lophsprc,
                     LLoutput,Lopusrgd,Lopusrhs,Lopgrdh)';
        if scalmiss(g);
            if not trapchk(4);
                errorlog "gradient calculation failed";
            endif;
            ret = g;
            goto A98;
        endif;
        isctu = 1;
    endif;
    if  iter >= Lopmiter;
        ret = error(2);
        goto A98;
    elseif ethsec(ttime,date)/6000 > Lopmtime;
        ret = error(11);
        goto A98;
    endif;

/*  test for convergence  */

    dfct = f0-vof;
    relgrad = (abs(g).*maxc(abs(x0)'|ones(1,rows(x0))))/maxc(abs(vof)|1);
    if abs(g) < smallval or relgrad < Lopgtol;
        ret = error(0);
        goto A98;
    endif;

    if Lopkey;
        gosub help;
    endif;

    if (midalg or midstp);
        if ((abs(100*(dfct/f0)) < Lopdfct) and iter > 1) or
            iter >= Lopditer or skip;
            chiter = 0;
            if midalg;
                k1 = Lopalgr;
                Lopalgr = midalg;
                midalg = 0;
                if Lopalgr >= 2 and Lopalgr <= 5;
                    h = eye(np)*maxc(sqrt(abs(vof))|1);
                endif;
            endif;
            if midstp;
                Lopstep = midstp;
                midstp = 0;
            endif;
            if midhs;
                gosub hssn;
                    pop h;
                isctu = 0;
                chiter = 0;
            endif;
        endif;
    endif;
    if isctu and midhs == 2 and chiter > Lopditer;
        gosub hssn;
            pop h;
        isctu = 0;
        chiter = 0;
    else;
        chiter = chiter + 1;
    endif;
    if LLoutput == 2;
        scroll w1;
    endif;
    goto A0;
A98:

    /* ******************** End of iteration loop ****************** */

    if Lopalgr == 1 or Lopalgr == 6;
        hsmeth = 0 $+ "NOHESS";
    elseif Lopalgr == 2 or Lopalgr == 3 or Lopalgr == 4;
        hsmeth = 0 $+ "SECANT";
    elseif Lopalgr == 5;
        hsmeth = 0 $+ "HESS";
    else;
        hsmeth = 0 $+ "NOHESS";
    endif;
    Lopitdta = (iter+1)|(ethsec(ttime,date)/6000)|hsmeth;
    if LLoutput == 2;
        cls;
    endif;
    goto OUT(x,vof,g,ret);

HELP:

    ky = key;
    do while ky;
    A5:
        if ky == 1030 or ky == 97 or ky == 65;      /* ALT A */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopalgr = ";
                locate 8,13;
                printdos ftos(Lopalgr,"%*.*lf",1,0);
                locate 10,4;
                printdos " = 1, steepest descent";
                locate 11,4;
                printdos " = 2, BFGS";
                locate 12,4;
                printdos " = 3, Modified BFGS";
                locate 13,4;
                printdos " = 4, DFP";
                locate 14,4;
                printdos " = 5, Newton-Raphson";
                locate 15,4;
                printdos " = 6, PRCG";
                locate 17,4;
                k1 = Lopalgr;
                printdos "Enter new value: ";
            else;
                print;
                print "Lopalgr = ";;
                print ftos(Lopalgr,"%*.*lf",1,0);
                print;
                print "    = 1, steepest descent";
                print "    = 2, BFGS";
                print "    = 3, Modified BFGS";
                print "    = 4, DFP";
                print "    = 5, Newton-Raphson";
                print "    = 6, PRCG";
                print;
                k1 = Lopalgr;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopalgr = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1120 or ky == 49;          /* ALT 1 */
            Lopalgr = 1;
            h = 1;
        elseif ky == 1121 or ky == 50;          /* ALT 2 */
            Lopalgr = 2;
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        elseif ky == 1122 or ky == 51;          /* ALT 3 */
            Lopalgr = 3;
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        elseif ky == 1123 or ky == 52;          /* ALT 4 */
            Lopalgr = 4;
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        elseif ky == 1124 or ky == 53;          /* ALT 5 */
            Lopalgr = 5;
            h = eye(np)*maxc(sqrt(abs(vof))|1);
        elseif ky == 1125 or ky == 54;          /* ALT 6 */
            Lopalgr = 6;
            h = 1;
        elseif ky == 1046 or upper(chrs(ky)) $== "C";          /* ALT C */
            ret = error(1);
            goto A98;       /* force convergence */
        elseif ky == 1018 or upper(chrs(ky)) $== "E";          /* ALT E */
            k1 = 1;
            if LLoutput == 2;
                scroll w5;
                locate 5,2;
                printdos "EDIT PARAMETER VECTOR";
                locate 5,28;
                printdos "<> <> to move   <ENTER> to select   <Q> to quit";
            A3_1:
                locate 6,4;
                printdos "par. no. ";
                printdos ftos(k1,"%*.*lf",3,0);
                printdos "  old value ";
                printdos ftos(x0[k1],"%*.*lf",10,6);
                ky = 0;
                do until ky == 81 or ky == 113;        /* Q or q */
                    ky = key;
                    if ky == 1072;
                        k1 = maxc(1|k1-1);
                        goto A3_1;
                    elseif ky == 1080;
                        k1 = minc(rows(x0)|k1+1);
                        goto A3_1;
                    elseif ky == 13;
                        locate 6,40;
                        printdos "new value ";
                        call csrtype(1);
                        x0[k1] = con(1,1);
                        call csrtype(0);
                        scroll w10;
                        goto A3_1;
                    endif;
                endo;
            else;
            A3_2:
                print;
                print " EDIT PARAMETER VECTOR    ";;
                print "<Bksp><Sp> to move   <ENTER> to select   <Q> to quit";
                print;
            A3_3:
                print "par. no. ";;
                print ftos(k1,"%*.*lf",3,0);;
                print "  old value ";;
                print ftos(x0[k1],"%*.*lf",10,6);;
                ky = 0;
                do until ky == 81 or ky == 113;        /* Q or q */
                    ky = key;
                    if ky == 8;                          /* Bksp */
                        k1 = maxc(1|k1-1);
                        print;
                        goto A3_3;
                    elseif ky == 32;                     /* Space */
                        k1 = minc(rows(x0)|k1+1);
                        print;
                        goto A3_3;
                    elseif ky == 13;                     /* Enter */
                        print "  new value ";
                        call csrtype(1);
                        x0[k1] = stof(cons);
                        print;
                        call csrtype(0);
                        goto A3_2;
                    endif;
                endo;
                print;
                print;
            endif;
            x = x0;
            gosub rebox;
            if LLoutput == 2;
                scroll w1;
            endif;
            if Lopalgr <= 4 or Lopalgr == 6;
                if LLoutput == 2;
                    locate 2,1;
                    printdos "H set to identity matrix";
                else;
                    print;
                    print "H set to identity matrix";
                    print;
                endif;
                if Lopalgr == 1 or Lopalgr == 6;
                    h = 1;
                else;
                    h = eye(np)*maxc(sqrt(abs(vof))|1);
                endif;

            elseif Lopalgr == 5;
                gosub hssn;
                    pop h;
            endif;
            if LLoutput == 2;
                scroll w1;
            endif;

        elseif ky == 1033 or upper(chrs(ky)) $== "F";          /* ALT F */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopdfct, change in function criterion";
                printdos ftos(Lopdfct,"%*.*lf",10,4);
                locate 10,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Lopdfct, change in function criterion";;
                print ftos(Lopdfct,"%*.*lf",10,4);
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopdfct = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1034 or upper(chrs(ky)) $== "G";          /* ALT G */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopgdmd = ";
                printdos ftos(Lopgdmd,"%*.*lf",1,0);
                locate 10,4;
                printdos " = 0, central difference method";
                locate 11,4;
                printdos " = 1, forward difference method";
                locate 12,4;
                printdos " = 2, forward difference, Richardson extrapolation"\
                    " method";
                locate 15,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Lopgdmd = ";;
                print ftos(Lopgdmd,"%*.*lf",1,0);
                print;
                print "    = 0, central difference method";
                print "    = 1, forward difference method";
                print "    = 2, forward difference, Richardson extrapolation"\
                         " method";
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopgdmd = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1023 or upper(chrs(ky)) $== "I";          /* ALT I */
            if LLoutput == 2;
                scroll w1;
            endif;
            if Lopalgr/=1;
                gosub hssn;
                    pop h;
            endif;
        elseif ky == 1050 or upper(chrs(ky)) $== "M";          /* ALT M */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Maximum number of backsteps = ";
                printdos ftos(Lopmbkst,"%*.*lf",1,0);
                locate 10,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Maximum number of backsteps = ";;
                print ftos(Lopmbkst,"%*.*lf",1,0);
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopmbkst = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1024 or upper(chrs(ky)) $== "O";          /* ALT O */
RETRY:
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "current setting, __OUTPUT = ";
                printdos ftos(LLoutput,"%*.*lf",3,0);
                locate 10,4;
                printdos "__OUTPUT = 0   no output";
                locate 11,4;
                printdos "__OUTPUT = 1   output suitable for file or printer";
                locate 12,4;
                printdos "__OUTPUT = 2   output suitable for screen only"\
                         " (ANSI.SYS required)";
                locate 14,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   current setting, __OUTPUT = ";;
                print ftos(LLoutput,"%*.*lf",3,0);
                print;
                print "__OUTPUT = 0   no output";
                print "__OUTPUT = 1   output suitable for file or printer";
                print "__OUTPUT = 2   output suitable for screen only"\
                      " (ANSI.SYS required)";
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                tmpop = stof(k0);
                if (tmpop /= 0 and tmpop /= 1 and tmpop /= 2) or
                    rows(tmpop) /= 1 or cols(tmpop) /= 1 or iscplx(tmpop);
                    if LLoutput == 2;
                        print "\007";
                    else;
                        print "Input error!";
                    endif;
                    goto retry;
                endif;
                if  tmpop == 2;             /* SPARC kludge */
                    print "2 not supported on SPARC yet, resetting to 1";
                    tmpop = 1;
                endif;
                LLoutput = tmpop;
            endif;
            gosub rebox;
            if LLoutput == 0;
                print "Output turned off.  Press Alt-O to reinstate.";
            endif;
            goto A9;
        elseif ky == 1025 or upper(chrs(ky)) $== "P";          /* ALT P */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopditer no. of iters for change in algorithm = ";
                printdos ftos(Lopditer,"%*.*lf",3,0);
                locate 10,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Lopditer no. of iters for change in algorithm = ";;
                print ftos(Lopditer,"%*.*lf",3,0);
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopditer = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1031 or upper(chrs(ky)) $== "S";          /* ALT S */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopstep = ";
                printdos ftos(Lopstep,"%*.*lf",1,0);
                locate 10,4;
                printdos " = 1, Step Length = 1";
                locate 11,4;
                printdos " = 2, STEPBT (Cubic, Quadratic)";
                locate 12,4;
                printdos " = 3, HALF    ";
                locate 13,4;
                printdos " = 4, BRENT   ";
                locate 15,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Lopstep = ";;
                print ftos(Lopstep,"%*.*lf",1,0);
                print;
                print "    = 1, Step Length = 1";
                print "    = 2, STEPBT (Cubic, Quadratic)";
                print "    = 3, HALF    ";
                print "    = 4, BRENT   ";
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopstep = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 33;    /* SHIFT 1 */
            Lopstep = 1;
        elseif ky == 64;    /* SHIFT 2 */
            Lopstep = 2;
        elseif ky == 35;    /* SHIFT 3 */
            Lopstep = 3;
        elseif ky == 36;    /* SHIFT 4 */
            Lopstep = 4;
        elseif ky == 1020 or upper(chrs(ky)) $== "T";          /* ALT T */
            skip = 1;
        elseif ky == 1047 or upper(chrs(ky)) $== "V";          /* ALT V */
            if LLoutput == 2;
                scroll w2;
                locate 8,4;
                printdos "Lopgtol gradient convergence criterion = ";
                printdos ftos(Lopgtol,"%*.*lf",10,6);
                locate 10,4;
                printdos "Enter new value: ";
            else;
                print;
                print "   Lopgtol gradient convergence criterion = ";;
                print ftos(Lopgtol,"%*.*lf",10,6);
                print;
                print "   Enter new value: ";;
            endif;
            call csrtype(1);
            k0 = cons;
            if LLoutput /= 2;
                print;
            endif;
            call csrtype(0);
            if k0 $/= "";
                Lopgtol = stof(k0);
            endif;
            gosub rebox;
        elseif ky == 1035 or upper(chrs(ky)) $== "H";          /* ALT H */
            call csrtype(0);
            if LLoutput == 2;
                scroll w5;
                locate 7,4;
                printdos "OPTIMIZATION SWITCHES";
                k0 = 196*ones(1,30);
                locate 8,4;
                printdos chrs(218~k0~194~k0~191);
                k1 = 0;
                do until k1 == 6;
                    k1 = k1+1;
                    locate 7+k1,4;
                    printdos chrs(179~(32*ones(1,30))~179~(32*ones(1,30))~179);
                endo;
                locate 8+k1,4;
                printdos chrs(192~k0~193~k0~217);
                locate 8,6;
                printdos "ALT I  Compute Hessian       ";
                locate 9,6;
                printdos "ALT G  Gradient Method       ";
                locate 10,6;
                printdos "ALT V  LopGTOL               ";
                locate 11,6;
                printdos "ALT F  LopDFCT               ";
                locate 12,6;
                printdos "ALT P  LopDITER              ";
                locate 13,6;
                printdos "ALT O  __OUTPUT              ";
                locate 8,37;
                printdos "ALT M  Maximum Backstep      ";
                locate 9,37;
                printdos "ALT E  Edit Parameter Vector ";
                locate 10,37;
                printdos "ALT C  Force Convergence     ";
                locate 11,37;
                printdos "ALT A  Algorithm             ";
                locate 12,37;
                printdos "ALT S  Step Length Method    ";
                locate 13,37;
                printdos "ALT T  Force Mid-method      ";
            else;
                print;
                print "OPTIMIZATION SWITCHES";
                print "---------------------";
                print "     ALT I  Compute Hessian       | ";;
                print "ALT M  Maximum Backstep      ";

                print "     ALT G  Gradient Method       | ";;
                print "ALT E  Edit Parameter Vector ";

                print "     ALT V  _MLGTOL               | ";;
                print "ALT C  Force Convergence     ";

                print "     ALT F  _MLDFCT               | ";;
                print "ALT A  Algorithm             ";

                print "     ALT P  _MLDITER              | ";;
                print "ALT S  Step Length Method    ";

                print "     ALT O  __OUTPUT              | ";;
                print "ALT T  Force Mid-method      ";

                print;
            endif;

            ky = key;
            do until ky;
                ky = key;
            endo;
            if LLoutput == 2 and
                ky /= 1047 and (not upper(chrs(ky)) $== "V") and
                ky /= 1050 and (not upper(chrs(ky)) $== "M") and
                ky /= 1020 and (not upper(chrs(ky)) $== "T") and
                ky /= 1034 and (not upper(chrs(ky)) $== "G") and
                ky /= 1025 and (not upper(chrs(ky)) $== "P") and
                ky /= 1033 and (not upper(chrs(ky)) $== "F");
                gosub parbox;
            endif;
            goto A5;
        elseif ky == 1073 or ky == 1081;    /* PgUp PgDn */
            if ky == 1081 and (np-pg*48) gt 0;
                pg = pg+1;
            elseif ky == 1073 and pg gt 1;
                pg = pg-1;
            else;
                goto A9;
            endif;
        endif;
    A9:

        ky = key;

    endo;

    return;

BAR:
    scroll w3;
    scroll w4;
    printdos "\27[7m";
    locate 1,1;
    printdos w0;
    locate 3,2;
    printdos "ITER      ";
    locate 3,19;
    printdos "FUNCTION: ";
    locate 4,2;
    printdos "TIME/ITER: ";
    locate 3,41;
    printdos "ALGORITHM: ";
    locate 3,60;
    printdos "STEP: ";
    locate 4,60;
    printdos "STEPSIZE: ";
    locate 3,52;
    printdos algrm[Lopalgr];
    locate 3,66;
    printdos stepm[Lopstep];
    locate 4,19;
    printdos "DF/ITER: ";
    locate 4,41;
    printdos "BACKSTEPS: ";
    locate 25,1;
    printdos "  ALT-H HELP  ";
    printdos "\27[0m";
    return;

PARBOX:

    scroll w5;
    locate 6,6;
    printdos "parameters/relative gradient";
    lf = (pg-1)*48+1;
    ll = minc(48|(np-(pg-1)*48));
    lr = ceil(ll/3);
    locate 7,1;
    printdos w6;
    k1 = 1;
    do until k1 gt lr;
        locate 7+k1,1;
        printdos w7;
        k1 = k1+1;
    endo;
    locate 7+k1,1;
    printdos w8;
    return;

REBOX:

    if LLoutput == 2;
        scroll w2;
        gosub parbox;
        locate 6,6;
        printdos "parameters/gradient";
    endif;
    return;

HSSN:

    h0 = {.};
    if LLoutput == 2;
        locate 2,15;
        printdos "hessian";
    endif;
    h0 = _deriv2(x0,2,&fnct,Lopgdmd,Lopgdprc,Lophsprc,LLoutput,
                  Lopusrgd,Lopusrhs,Lopgrdh);
    if scalmiss(h0);
        if not trapchk(4);
            if LLoutput == 2;
                locate 2,1;
            endif;
            errorlog "Hessian calculation failed";
        endif;
        h0 = eye(np)*maxc(sqrt(abs(vof))|1);
    endif;
    return(h0);

CHK1:
        pop c1;
        pop c0;
    if c0 and not scalmiss(c0);
       if c0 == 7;
          c0 = 5;
       endif;
       c0 = maxc(c0|1);
       c0 = minc(c0|6);
    else;
       c0 = c1;
    endif;
    return(c0);

CHK2:
        pop c1;
        pop c0;
    if c0 and not scalmiss(c0);
      c0 = maxc(c0|1);
      if c0 == 5 or c0 == 7;
        c0 = 1;
      endif;
      if c0 == 7;
        c0 = 4;
      endif;
    else;
      c0 = c1;
    endif;
    return(c0);

CHK3:
        pop c1;
        pop c0;

    if c0 and not scalmiss(c0);
        if c0 > 0;
           c0 = 1;
        else;
           c0 = 0;
        endif;
    else;
       c0 = c1;
    endif;
    return(c0);

CHK4:
        pop c0;
    if rows(c0) == 1 and cols(c0) == 1;
       if c0 > 1;
          c0 = 1;
       elseif c0 <= 1;
          c0 = 0;
       endif;
    endif;
    return(c0);

OUT:
     pop ret;
     pop g;
     pop vof;
     pop x;
     if not scalmiss(h);
         if Lopalgr < 5 and Lopalgr > 1;
            Lopfhess = h'h;
         elseif Lopalgr == 5;
            Lopfhess = h;
         else;
            Lopfhess = error(0);
         endif;
     endif;
     call ndpcntrl(old,0xffff);
     ndpclex;
retp(x,vof,g,scalerr(ret),Lopfhess,Lopitdta);
endp;

proc _strsch(a,b);
    local old;
    old = ndpcntrl(0,0);
    call ndpcntrl(0x0002,0x0002);
    a = packr(indcv(missrv(miss(stof(a),1),"ONE"),b));
    ndpclex;
    call ndpcntrl(old,0xffff);
    retp(a);
endp;

/*-----------------------------------------------------*/
/*   PROC SCTU                                         */
/* This computes G & updates to inverse Hessian        */
/*-----------------------------------------------------*/

proc(3) = _sctu2(x,vof,smallval,g,h,dx,d,s,&fnct,Lopalgr,LLoutput,
                Lopgdmd,Lopgdprc,Lophsprc,Lopusrgd,Lopusrhs,Lopgrdh);

/*------- LOCALS --------*/
    local v1, v2, g0, h1, w1, w2, beta;
    local fnct:proc;
    clear beta;
    h1 = h;

    /* --- Gradient at x --- */
    g0 = g;
    if LLoutput == 2;
        locate 2,15;
        printdos " Gradient";
    endif;
    g = _deriv2(x,1,&fnct,Lopgdmd,Lopgdprc,Lophsprc,LLoutput,
                 Lopusrgd,Lopusrhs,Lopgrdh)';
    if scalmiss(g);
        if not trapchk(4);
            errorlog "gradient calculation failed";
        endif;
        retp(g,h1,beta);
    endif;

    if abs(g) < smallval;
        retp(g,h1,beta);
    endif;

    if Lopalgr == 1;         /*  Steepest Descent  */

        retp(g,0,0);

    elseif Lopalgr == 2;        /* BFGS update */

        v1 = g'dx - g0'dx;
        if (v1 < 1e-22);
            h1 = eye(rows(x))*maxc(sqrt(abs(vof))|1);
        else;
            h1 = cholup(h,(g-g0)/sqrt(v1));
            h1 = choldn(h1,g0/sqrt(-g0'd));
        endif;

    elseif Lopalgr == 3;    /* scale-free BFGS update */

        v1 = g'dx - g0'dx;
        if (v1 < 1e-222);
            h1 = eye(rows(x))*maxc(sqrt(abs(vof))|1);
        else;
            w1 = sqrt(-v1/(g0'd))/s;
            h1 = cholup(w1*h,(g-g0)/sqrt(v1));
            h1 = choldn(h1,g0*w1);
        endif;

    elseif Lopalgr == 4;    /* DFP update */

        v1 = g'dx - g0'dx;
        if (v1 < 1e-222);
            h1 = eye(rows(x))*maxc(sqrt(abs(vof))|1);
        else;
            v2 = sqrt(-g0'd);
            w1 = sqrt(v1);
            h1 = cholup(h,g/w1-g0/w1);
            w2 = s*(v2/v1);
            h1 = cholup(h1,w2*g-w2*g0+g0/v2);
            h1 = choldn(h1,g0/v2);
        endif;

    elseif Lopalgr == 5;      /* NEWTON-RAPHSON if 5  */

        if LLoutput == 2;
            locate 2,15;
            printdos " Hessian";
        endif;
        h1 = _deriv2(x,2,&fnct,Lopgdmd,Lopgdprc,Lophsprc,LLoutput,
                    Lopusrgd,Lopusrhs,Lopgrdh)';
        if scalmiss(h1);
            if not trapchk(4);
                if LLoutput == 2;
                    locate 2,15;
                endif;
                errorlog "Hessian calculation failed";
            endif;
            h1 = eye(rows(x))*maxc(sqrt(abs(vof))|1);
        endif;
    elseif Lopalgr == 6;
        h1 = 1;
        v1 = g0'*g0;
        beta = (g'g/v1)*dx - (g0'g/v1)*d;
    endif;
    retp(g,h1,beta);
endp;

proc(2) = _stepl2(g,vof,x0,d,&fnct,Lopstep,Lopusrch,Lopmbkst,LLoutput,Lopmxtry);
    local s, rs, ret, bksteps;
    local x,y,w;
    local fnct:proc;
    let w = 2 1 2 9 0 7;

    clear ret;
    bksteps = -1;
    if Lopstep == 2;
        gosub fct(x0+d); pop rs;
        { s,ret,bksteps } = _stepb2(rs,g,vof,x0,d,&fnct,Lopmbkst,LLoutput);
    elseif Lopstep == 3;
        { s,ret,bksteps } = _half2(vof,x0,d,&fnct,Lopmxtry,LLoutput);
    elseif Lopstep == 4;
        { s,ret,bksteps } = _brent2(vof,x0,d,1e-2,&fnct,Lopmxtry,LLoutput);
    else;
        s = 1;
        gosub fct(x0+d); pop rs;
        if rs > vof;
            ret = 1;
        endif;
    endif;

    if ret == 1 and Lopstep /= 4;    /* not successful */
        if LLoutput == 2;
            locate 3,66;
            printdos ("\27[7m" $+ "BRENT " $+ "\27[0m");
        endif;
        { s,ret,bksteps } = _brent2(vof,x0,d,1e-2,&fnct,Lopmxtry,LLoutput);
    endif;

    if ret == 1 and Lopstep /= 3;    /* still not successful */
        if LLoutput == 2;
            locate 3,66;
            printdos ("\27[7m" $+ "HALF  " $+ "\27[0m");
        endif;
        { s,ret,bksteps } = _half2(vof,x0,d,&fnct,Lopmxtry,LLoutput);
    endif;

    if ret == 1 and Lopusrch;
        { s,ret } = _usrsch2(s,x0,d,vof,&fnct,LLoutput);
    endif;

    if ret == 1;
        s = error(6);
    endif;

    retp(s,bksteps);

FCT:
    pop x;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    y = fnct(x);
    if LLoutput == 2;
        scroll w;
    endif;
    if scalmiss(y) or
        y == __INFp or
        y == __INFn or
        y == __INDEFn or
        y == __INDEFp;
        retp(error(3),bksteps);
    endif;
return(y);

endp;



proc _deriv2(x,ind,&fct,gdmd,gdprc,hsprc,LLoutput,Lopusrgd,Lopusrhs,Lopgrdh);
    local a,i,gdproc,hsproc;
    local fct:proc, Lopusrgd:proc;

    if LLoutput == 2;
        locate 2,15;
        if ind == 1;
            printdos " Gradient";
        elseif ind == 2;
            printdos " Hessian";
        endif;
    endif;

    if hsprc /= 0;
        hsproc = hsprc;
        local hsproc:proc;
    endif;

    if gdprc /= 0;
        gdproc = gdprc;
        local gdproc:proc;
    endif;
    if ind == 1;
        if gdprc /= 0;
            a = gdproc(x);
        else;
            if gdmd == 0;
                a = _gdcd(&fct,x,Lopgrdh);
            elseif gdmd == 1;
                a = _gdfd(&fct,x,Lopgrdh);
            else;
                a = Lopusrgd(&fct,x);
            endif;
        endif;
    elseif ind == 2;
        if hsprc /= 0;
            a = hsproc(x);
        else;
            if gdprc /= 0;
                if gdmd == 0;
                    a = _gdcd(gdprc,x,Lopgrdh);
                elseif gdmd == 1;
                    a = _gdfd(gdprc,x,Lopgrdh);
                else;
                    a = Lopusrgd(gdprc,x);
                endif;
            else;
                if Lopusrhs /= 0;
                    local Lopusrhs:proc;
                    a = Lopusrhs(&fct,x);
                else;
                    a = _hssp(&fct,x,Lopgrdh);
                endif;
            endif;
        endif;
    endif;
    if scalmiss(a);
        if not trapchk(4);
            if LLoutput == 2;
                locate 2,15;
            endif;
            errorlog "Calculation of derivatives failed";
        endif;
        if ind == 1;
           retp(error(4));
        else;
           retp(error(5));
        endif;
    endif;
    if LLoutput == 2;
        let i = 2 15 2 25 0 7;
        scroll i;
    endif;
    retp(a);
endp;




proc _gdfd(f,x0,Lopgrdh);
    local f:proc;
    local n, k, grdd, dh, ax0, xdh, arg, dax0, i, f0, eps;

    f0 = f(x0);
    n = rows(f0);
    k = rows(x0);
    grdd = zeros(n,k);
    eps = 6.0554544523933429e-6;

/* Computation of stepsize (dh) for gradient */

    if Lopgrdh /= 0;
        dh = Lopgrdh;
    else;
        ax0 = abs(x0);
        if x0 /= 0;
            dax0 = x0 ./ ax0;
        else;
            dax0 = 1;
        endif;
        dh = eps*maxc((ax0~(1e-2)*ones(k,1))').*dax0;
    endif;

    xdh = x0+dh;
    dh = xdh-x0;    /* This increases precision slightly */
    arg = diagrv(reshape(x0,k,k)',xdh);

    i = 1;
    do until i > k;
        grdd[.,i] = f(submat(arg,0,i))';
        i = i+1;
    endo;
    grdd = (grdd-f0)./(dh');
    retp(grdd);
endp;


proc _gdcd(f,x,Lopgrdh);
    local f:proc,k,ax0,dax0,dh,xdh,argplus,argminus,i,grdd,eps;
    k = rows(x);
    eps = 6.0554544523933429e-6;

    if Lopgrdh /= 0;
        dh = Lopgrdh;
    else;
        ax0 = abs(x);
        if x /= 0;
            dax0 = x./ax0;
        else;
            dax0 = 1;
        endif;
        dh = eps*maxc((ax0~(1e-2)*ones(k,1))').*dax0;
    endif;

    xdh = x+dh;
    dh = xdh-x;     /* This increases precision slightly */
    argplus = diagrv(reshape(x,k,k)',xdh);
    argminus = diagrv(reshape(x,k,k)',x-dh);

    i = 0;
    do until i == k;
        i = i + 1;
        if i == 1;
            grdd = f(argplus[.,1])'-f(argminus[.,1])';
        else;
            grdd = grdd~(f(argplus[.,i])'-f(argminus[.,i])');
        endif;
    endo;
    retp(grdd./(2*dh'));
endp;



proc _hssp(&f,x0,Lopgrdh);
    local k, hessian, grdd, ax0, dax0, dh, xdh, ee, f0, i, j, eps;
    local f:proc;

    /* check for complex input */
    if iscplx(x0);
        if hasimag(x0);
            errorlog "ERROR: Not implemented for complex matrices.";
            end;
        else;
            x0 = real(x0);
        endif;
    endif;

/* initializations */
    k = rows(x0);
    hessian = zeros(k,k);
    grdd = zeros(k,1);
    eps = 6.0554544523933429e-6;

/* Computation of stepsize (dh) */
    if Lopgrdh /= 0;
        dh = Lopgrdh;
    else;
        ax0 = abs(x0);
        if x0 /= 0;
            dax0 = x0./ax0;
        else;
            dax0 = 1;
        endif;
        dh = eps*maxc((ax0~(1e-2)*ones(k,1))').*dax0;
    endif;

    xdh = x0+dh;
    dh = xdh-x0;    /* This increases precision slightly */
    ee = eye(k).*dh;

/* Computation of f0=f(x0) */
    f0 = f(x0);

/* Compute forward step */
    i = 1;
    do until i > k;

        grdd[i,1] = f(x0+ee[.,i]);

        i = i+1;
    endo;

/* Compute "double" forward step */
    i = 1;
    do until i > k;
        j = i;
        do until j > k;

            hessian[i,j] = f(x0+(ee[.,i]+ee[.,j]));
            if i /= j;
                hessian[j,i] = hessian[i,j];
            endif;

            j = j+1;
        endo;
        i = i+1;
    endo;

    retp( ( ( (hessian - grdd) - grdd') + f0) ./ (dh.*dh') );
endp;



proc(3) = _stepb2(r1,g,vof,x0,d,&fnct,mbkst,LLoutput);
    local delta,ub,lb, ret, i, cdelta, dg, s, g1, r2, rs, sprev, s2prev,
        tt, rprev, r2prev, sprev2, s2prev2, sp2, dsprev, vv, zz, ab, a, b,
        qv;

    local x,y,w;
    local fnct:proc;
    let w = 2 1 2 9 0 7;

/* --------------------- Initializations -------------------------  */
    delta = 1e-4;           /* This can be changed, and doing so may help  */
            /* speed convergence -- it must remain within the interval  */
            /* (0,1/2) */
    ub = 0.5;       /* Upper bound on acceptable reduction in s. */
    lb = 0.1;       /* Lower bound on acceptable reduction in s. */

    ret = 1;        /* If 0, then satisfactory value found; else 1.  */
    i = 0;          /* This counts # of backsteps taken. */

    cdelta = 1-delta;

    dg = d'*g;

        /* ------------------- Try s=1 -------------------------- */
    s = 1;
    tt = s*dg;
    g1 = r1/tt-vof/tt;
    if g1>=delta;
        if r1 > vof;
            retp(s,1,0);
        else;
            retp(s,0,0);
        endif;
    endif;
    i = 1;
    s = -dg/(2*(r1-vof-dg));
    s = maxc(s|lb);
    gosub fct(x0+s'.*d); pop r2;
    tt = s*dg;
    g1 = r2/tt-vof/tt;

    if g1>=delta and g1<=cdelta;
        if r2 > vof;
            retp(s,1,1);
        else;
            retp(s,0,1);
        endif;
    endif;
    sprev = s;
    s2prev = 1;
    rprev = r2;
    r2prev = r1;
    rs = r2;

    i = 2;
    do until i == mbkst;

        sprev2 = sprev*sprev;
        s2prev2 = s2prev*s2prev;
        sp2 = sprev2~s2prev2;
        dsprev = sprev-s2prev;

        vv = (1~-1|-s2prev~sprev);
        vv = vv./sp2;
        zz = (rprev-vof-dg*sprev)|(r2prev-vof-dg*s2prev);
        ab = (1/dsprev)*vv*zz;
        a = ab[1,1];
        b = ab[2,1];

        if a == 0;          /* Cubic is actually a Quadratic in this case. */
            s = -dg/(2*b);
        else;
            qv = b*b - 3*a*dg;
            if qv < 0;
                if rs > vof;
                    retp(s,1,i);
                else;
                    retp(s,ret,i);
                endif;
            endif;          /* terminate if not real root */
            tt = 3*a;
            s = -b/tt + sqrt(qv)/tt;
        endif;

        if s > ub*sprev;
            s = ub*sprev;
        elseif s < lb*sprev;
            s = lb*sprev;
        endif;

        gosub fct(x0+s'.*d); pop rs;
        tt = s*dg;
        g1 = rs/tt-vof/tt;
        if g1 >= delta and g1 <= cdelta;
            if rs > vof;
                retp(s,1,i);
            else;
                retp(s,0,i);
            endif;
        endif;

        s2prev = sprev;
        sprev = s;
        r2prev = rprev;
        rprev = rs;
        i = i+1;
    endo;
    if rs > vof;
        retp(s,1,i);
    else;
        retp(s,ret,i);
    endif;

FCT:
    pop x;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    y = fnct(x);
    if LLoutput == 2;
        scroll w;
    endif;
    if scalmiss(y) or
        y == __INFp or
        y == __INFn or
        y == __INDEFn or
        y == __INDEFp;
        retp(error(3),1,i);
    endif;
return(y);


endp;

proc(2) = _usrsch2(s,x0,d,vof,&fnct,LLoutput);
    local rs,w,w1,w2,eps,ky;
    local fnct:proc;
    let w = 2 30 2 80 0 7;
    let w1 = 5 1 24 80 0 7;
    let w2 = 2 1 2 9 0 7;
    if LLoutput == 2;
        scroll w1;
        locate 2,30;
        printdos "entering user search";
        locate 6,4;
        printdos "initial function value ";
        printdos ftos(vof,"%*.*lf",16,8);
    else;
        print;
        print "entering user search";
        print;
        print "initial function value ";;
        print ftos(vof,"%*.*lf",16,8);
        print;
    endif;
    eps = .1*s;
    clear ky;
A0:

    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    else;
        print " Function";
        print;
    endif;
    rs = fnct(x0+s*d);
    if LLoutput == 2;
        scroll w2;
    endif;
    if scalmiss(rs);
      retp(error(3),1);
    endif;
    if LLoutput == 2;
        locate 8,4;
        printdos "stepsize ";
        printdos ftos(s,"%*.*lf",16,8);
        printdos "    eps = ";
        printdos ftos(eps,"%*.*lf",16,8);
        locate 10,4;
        printdos "function ";
        printdos ftos(rs,"%*.*lf",16,8);
    else;
        print "   stepsize ";;
        print ftos(s,"%*.*lf",16,8);;
        print "    eps = ";;
        print ftos(eps,"%*.*lf",16,8);
        print;
        print "   function ";;
        print ftos(rs,"%*.*lf",16,8);
        print;
    endif;
    ky = key;
    do until ky;
        ky = key;
    endo;
    if ky == 27;
        if LLoutput == 2;
            scroll w;
            scroll w1;
        endif;
        retp(s,0);
    elseif ky == 1072;
        s = s+eps;
    elseif ky == 1080;
        s = s-eps;
    elseif ky == 56;
        eps = eps*10;
    elseif ky == 50;
        eps = eps/10;
    endif;
    goto A0;
endp;


proc(2) = _brackt2(f0,x,d,&fnct,mxtry,LLoutput);

    local g,r,l1,f1,l2,f2,try;
    local t,y,w,w1;
    local fnct:proc;
    let w = 2 1 2 9 0 7;
    let w1 = 2 15 2 26 0 7;
    if LLoutput == 2;
        locate 2,15;
        printdos "bracketing";
    endif;
    g = 0.5*sqrt(5) - 0.5;
    r = 2 - g;
    l1 = g;
    gosub fct(x+l1*d); pop f1;
    try = 0;
    do until try > mxtry;
       try = try + 1;
       if f1 > f0;
           l2 = g*l1;
           gosub fct(x+l2*d); pop f2;
           if (f2 < f0);
               if LLoutput == 2;
                    scroll w1;
               endif;
               retp(0,l1);
           else;
               l1 = l2;
               f1 = f2;
           endif;
       else;
           l1 = r*l1;
           gosub fct(x+l1*d); pop f1;
       endif;
    endo;
    retp(error(0),error(0));

FCT:
    pop t;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    y = fnct(t);
    if LLoutput == 2;
        scroll w;
    endif;
    if scalmiss(y) or
        y == __INFp or
        y == __INFn or
        y == __INDEFn or
        y == __INDEFp;
        retp(error(0),error(0));
    endif;
return(y);

endp;


proc(3) = _brent2(vof,x0,d0,tol,&fnct,mxtry,LLoutput);

    local ax,cx,iter;
    local a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w,fu,fv,fw,fx,x,tol3;
    local f,y,w1;
    local fnct:proc;
    let w1 = 2 1 2 9 0 7;

    eps = sqrt(2.2204460492503131e-16);
    c = 0.5*(3 - sqrt(5));

    {ax,cx} = _brackt2(vof,x0,d0,&fnct,mxtry,LLoutput);
    if scalmiss(ax);
        retp(ax,1,0);
    endif;
    a = minc(ax|cx);
    b = maxc(ax|cx);
    v = a + c*(b - a);
    w = v;
    x = v;
    e = 0;
    fx = vof;
    fv = fx;
    fw = fx;
    tol3 = tol/3;

    iter = 0;
    do until iter == mxtry;
        iter = iter + 1;
        xm = 0.5*(a + b);
        tol1 = eps*abs(x) + tol3;
        t2 = 2*tol1;
        if abs(x - xm) <= (t2 - .5*(b - a));
            retp(x,0,iter);
        endif;
        clear p,q,r;
        if abs(e) <= tol1;
            goto A40;
        endif;
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r;
        q = 2*(q-r);
        if q <= 0;
            q = -q;
        else;
            p = -p;
        endif;
        r = e;
        e = d;
        if (abs(p) >= abs(.5*q*r)); goto A40; endif;
        if (p <= q*(a-x)); goto A40; endif;
        if (p >= q*(b-x)); goto A40; endif;
        d = p/q;
        u = x + d;
        if ((u-a) < t2);
            if (xm - x) > 0;
                d = tol1;
            else;
                d = - tol1;
            endif;
        endif;
        if ((b-u) < t2);
            if (xm - x) > 0;
                d = tol1;
            else;
                d = - tol1;
            endif;
        endif;
        goto A50;
A40:
        if x < xm;
            e = b-x;
        else;
            e = a-x;
        endif;
        d = c*e;
A50:
        if abs(d) >= tol1;
            u = x + d;
        elseif d >= 0;
            u = x + tol1;
        else;
            u = x - tol1;
        endif;
        gosub fct(x0 + u*d0); pop fu;
        if fu <= fx;
            if u >= x;
                a = x;
            else;
                b = x;
            endif;
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        else;
            if u >= x;
                b = u;
            else;
                a = u;
            endif;
            if (fu <= fw) or (w == x);
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            elseif (fu <= fv) or (v == x) or (v == w);
                v = u;
                fv = fu;
            endif;
        endif;
    endo;
    retp(error(6),1,iter);

FCT:
    pop f;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    y = fnct(f);
    if LLoutput == 2;
        scroll w1;
    endif;
    if scalmiss(y) or
        y == __INFp or
        y == __INFn or
        y == __INDEFn or
        y == __INDEFp;
        retp(error(3),1,iter);
    endif;
return(y);

endp;


proc(3) = _half2(f,x,d,&fnct,mxtry,LLoutput);
   local ax,cx,bksteps,f1,t,y,w;
   local fnct:proc;
   let w = 2 1 2 9 0 7;

   {ax,cx} = _brackt2(f,x,d,&fnct,mxtry,LLoutput);
   bksteps = 0;
   do until bksteps > mxtry;
      bksteps = bksteps + 1;
      cx = ax + .5*(cx - ax);
      gosub fct(x+cx*d); pop f1;
      if f1 < f;
         retp(cx,0,bksteps);
      endif;
   endo;
   retp(miss(0,0),1,bksteps);

FCT:
    pop t;
    if LLoutput == 2;
        locate 2,1;
        printdos " Function";
    endif;
    y = fnct(t);
    if LLoutput == 2;
        scroll w;
    endif;
    if scalmiss(y) or
        y == __INFp or
        y == __INFn or
        y == __INDEFn or
        y == __INDEFp;
        retp(error(3),1,bksteps);
    endif;
return(y);

endp;
