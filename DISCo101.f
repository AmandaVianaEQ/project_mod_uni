C-------------------------------------------------------------------
       SUBROUTINE BANDIDX(INDEX,NEQ,NE,MU,ML,LINDEX,LVALUE)
       
        IMPLICIT NONE
        INTEGER NEQ(*),INDEX(*),I,J,NE,LINDEX,LVALUE,L,ML,MU
        INTEGER LENP
      
        LENP=(ML+MU+1)*NEQ(1)
        NE=(ML+MU+1)*NEQ(1)
c       LINDEX=3*NE+36*NEQ(1)+1+100
c       LVALUE=4*NE
        LINDEX=2*NE
        LVALUE=NE
        CALL NULLIDX(INDEX,LINDEX)
        L=0 
        DO I=1,NEQ(1)
          DO J=1,ML+MU+1
            L=L+1
            INDEX(L)=J+I-MU-1
            INDEX(L+LENP)=I
          END DO
        END DO
        RETURN
        END
      SUBROUTINE CFODE(METH,ELCO,TESCO)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------
CLLL. OPTIMIZE

      IMPLICIT NONE

      INTEGER METH,I,IB,NQ,NQM1,NQP1
      DOUBLE PRECISION ELCO(13,12),TESCO(3,12),AGAMQ,FNQ,FNQM1,PC(12),
     1                 PINT,RAGQ,RQFAC,RQ1FAC,TSIGN,XPIN
C---------------------------------------------------------------------
C CFODE IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
C NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS 
C GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.  
C THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH = 2.
C (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
C CFODE IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM, 
C AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
C                                                        
C THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.    
C THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF 
C ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A GENETRATING
C POLYNOMIAL, I.E.,                                
C     L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.  
C FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY     
C     DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) = 0.
C FOR THE BDF METHODS, L(X) IS GIVEN BY    
C     L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,       
C WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
C                                                   
C THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE  
C LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
C AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP   
C SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
C NQ + 1 IF K = 3.       
C---------------------------------------------------------------------

      GO TO (100,200),METH

 100  ELCO(1,1)=1.0D0
      ELCO(2,1)=1.0D0
      TESCO(1,1)=0.0D0
      TESCO(2,1)=2.0D0
      TESCO(1,2)=1.0D0
      TESCO(3,12)=0.0D0
      PC(1)=1.0D0
      RQFAC=1.0D0
      DO NQ=2,12

C---------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
C     P(X) = (X+1)*(X+2)*...*(X+NQ-1).      
C INITIALLY, P(X) = 1.   
C---------------------------------------------------------------------

        RQ1FAC=RQFAC
        RQFAC=RQFAC/DFLOAT(NQ)
        NQM1=NQ-1
        FNQM1=DFLOAT(NQM1)
        NQP1=NQ+1

C FORM COEFFICIENTS OF P(X)*(X+NQ-1). --------------------------------

        PC(NQ)=0.0D0
        DO IB=1,NQM1
          I=NQP1-IB
          PC(I)=PC(I-1)+FNQM1*PC(I)
        END DO
        PC(1)=FNQM1*PC(1)

C COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X). ---------------------

        PINT=PC(1)
        XPIN=PC(1)/2.0D0
        TSIGN=1.0D0
        DO I=2,NQ
          TSIGN=-TSIGN
          PINT=PINT+TSIGN*PC(I)/DFLOAT(I)
          XPIN=XPIN+TSIGN*PC(I)/DFLOAT(I+1)
        END DO

C STORE COEFFICIENTS IN ELCO AND TESCO. ------------------------------

        ELCO(1,NQ)=PINT*RQ1FAC
        ELCO(2,NQ)=1.0D0
        DO I=2,NQ
          ELCO(I+1,NQ)=RQ1FAC*PC(I)/DFLOAT(I)
        END DO
        AGAMQ=RQFAC*XPIN
        RAGQ=1.0D0/AGAMQ
        TESCO(2,NQ)=RAGQ
        IF (NQ.LT.12) TESCO(1,NQP1)=RAGQ*RQFAC/DFLOAT(NQP1)
        TESCO(3,NQM1)=RAGQ
      END DO
      RETURN
C    
 200  PC(1)=1.0D0
      RQ1FAC=1.0D0
      DO NQ=1,5

C---------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
C     P(X) = (X+1)*(X+2)*...*(X+NQ).
C INITIALLY, P(X) = 1.
C---------------------------------------------------------------------

        FNQ=DFLOAT(NQ)
        NQP1=NQ+1

C FORM COEFFICIENTS OF P(X)*(X+NQ). ----------------------------------

        PC(NQP1)=0.0D0
        DO IB=1,NQ
          I=NQ+2-IB
          PC(I)=PC(I-1)+FNQ*PC(I)
        END DO
        PC(1)=FNQ*PC(1)

C STORE COEFFICIENTS IN ELCO AND TESCO. ------------------------------

        DO I=1,NQP1
          ELCO(I,NQ)=PC(I)/PC(2)
        END DO
        ELCO(2,NQ)=1.0D0
        TESCO(1,NQ)=RQ1FAC
        TESCO(2,NQ)=DFLOAT(NQP1)/ELCO(1,NQ)
        TESCO(3,NQ)=DFLOAT(NQ+2)/ELCO(1,NQ)
        RQ1FAC=RQ1FAC/FNQ
      END DO
      RETURN
C----------------------- Fin de la routine CFODE ---------------------
      END
      DOUBLE PRECISION FUNCTION D1MACH()

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE IN DOUBLE
C PRECISION.  THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
C U SUCH THAT  1.0D0 + U .NE. 1.0D0 (IN DOUBLE PRECISION).
C---------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION U,COMP

      U=1.0D0
      COMP=1.0D0+U
      DO WHILE(COMP.NE.1.0D0)
        U=U*0.5D0
        COMP=1.0D0+U
      END DO
      D1MACH=U*2.0D0
      RETURN
C----------------------- Fin de la fonction D1MACH ---------------------        
      END
C---------------------------------------------------------------------
      SUBROUTINE DAXEY(N,DA,DX,INCX,DY,INCY)
C                                   
      IMPLICIT NONE
      DOUBLE PRECISION DX(*),DY(*),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C                    
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C                                                      
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1 
C          
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C                                          
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C                   
C                    
C        CLEAN-UP LOOP
C  
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DA*DX(I)
        DY(I + 1) =  DA*DX(I + 1)
        DY(I + 2) =  DA*DX(I + 2)
        DY(I + 3) =  DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C-------------------------------------------------------------------
       SUBROUTINE DENSIDX(INDEX,NEQ,NE,LINDEX,LVALUE)
       
        IMPLICIT NONE
        INTEGER NEQ(*),INDEX(*),I,J,NE,LINDEX,LVALUE

        NE=NEQ(1)*NEQ(1)
c       LINDEX=3*NE+36*NEQ(1)+1+100
c       LVALUE=4*NE
        LINDEX=2*NE
        LVALUE=NE
        CALL NULLIDX(INDEX,LINDEX)
        DO I=1,NEQ(1)
         DO J=1,NEQ(1)
          INDEX((I-1)*NEQ(1)+J)=J
          INDEX(NEQ(1)*NEQ(1)+(I-1)*NEQ(1)+J)=I
         END DO
        END DO
        RETURN
        END
      SUBROUTINE DISCO(RES,JAC,ADDA,NEQ,Y,YDOTI,T,TOUT,ITOL,RTOL,
     1                 ATOL,ITASK,ISTATE,IOPT,RW,LRW,IW,LIW,MF,
     2                 GEX,NG,JROOT,ALPHA,RPAR,LRP,IPAR,LIP)
C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL RES,JAC,ADDA,GEX

      INTEGER NEQ(*),ITOL,ITASK,ISTATE,IOPT,LRW,LIW,IW(*),MF,NG,
     1        JROOT(*),INDPRP,I,I1,I2,IER,IFLAG,IMXER,IRES,KGO,
     2        LENIW,LENRW,LENWM,LP,LYD0,ML,MORD(2),MU,MXHNL0,
     3        MXSTP0,IRFP,IRT,LENYH,LYHNEW,IPAR(*),LRP,LIP

      DOUBLE PRECISION CON,RPAR(*),HINI,T,TOUT,Y(*),YDOTI(*),RTOL(*),
     1        ATOL(*),RW(*),ATOLI,AYI,BIG,
     2        EWTI,H0,HMAX,HMX,RH,RTOLI,TCRIT,TDIST,TNEXT,TOL,TOLSF,TP,
     3        SIZE,SUM,W0,D1MACH,VNORM,ALPHA(*)

      LOGICAL IHIT


C-----------------------------------------------------------------------
C THE FOLLOWING INTERNAL COMMON BLOCK CONTAINS
C (A) VARIABLES WHICH ARE LOCAL TO ANY SUBROUTINE BUT WHOSE VALUES MUST
C     BE PRESERVED BETWEEN CALLS TO THE ROUTINE (OWN VARIABLES), AND
C (B) VARIABLES WHICH ARE COMMUNICATED BETWEEN SUBROUTINES.
C COMMON BLOCK LS0001 IS SHARED BY THE DISCo AND LSODE PACKAGES.
C THE STRUCTURE OF LS0001 IS AS FOLLOWS..  ALL REAL VARIABLES ARE
C LISTED FIRST, FOLLOWED BY ALL INTEGERS.  WITHIN EACH TYPE, THE
C VARIABLES ARE GROUPED WITH THOSE LOCAL TO SUBROUTINE DISCo FIRST,
C THEN THOSE LOCAL TO SUBROUTINE STODI, AND FINALLY THOSE USED
C FOR COMMUNICATION.  THE BLOCK IS DECLARED IN SUBROUTINES
C DISCo, INTDY, STODI, PREPJI, AND SOLSY.  GROUPS OF VARIABLES ARE
C REPLACED BY DUMMY ARRAYS IN THE COMMON DECLARATIONS IN ROUTINES
C WHERE THOSE VARIABLES ARE NOT USED.
C-----------------------------------------------------------------------

C------------------------------------------------------------------------
      INTEGER LTRET,LCONIT,LCRATE,LEL,LELCO,LHOLD,LRMAX,LTESCO,LCCMAX,
     1        LEL0,LH,LHMIN,LHMXI,LHU,LRC,LTN,LUROUND,LALPHA,LX2,LT0,
     2        LTLAST,LTOUTC,LCNTL,LRINFO,LMBAND,LNE,LCHOUT,LINDEX,
     3        LVALUE,LILLIN,LINIT,LLYH,LLEWT,LLACOR,LLSAVR,LLWM,LLIWM,
     4        LMXSTEP,LMXHNIL,LNHNIL,LNTREP,LNSLAST,LNYH,LIALTH,LIPUP,
     5        LLMAX,LMEO,LNQNYH,LNSLP,LICF,LIERPJ,LIERSL,LJCUR,LJSTART,
     6        LKFLAG,LL,LMETH,LMITER,LMAXORD,LMAXCOR,LMSBP,LMXNCF,LN,
     7        LNQ,LNST,LNRE,LNJE,LNQU,LLG0,LLG1,LLGX,LIMAX,LLAST,
     8        LIRFND,LITASKC,LNGC,LNGE,LMESFLG,LLUNIT,LKEEP,LINFO,
     9        LICNTL,LIDX

      DATA MORD(1),MORD(2)/12,5/,MXSTP0/500/,MXHNL0/10/

      
      LTRET=21
      LCONIT=LTRET+1
      LCRATE=LCONIT+1
      LEL=LCRATE+1
      LELCO=LEL+13
      LHOLD=LELCO+156
      LRMAX=LHOLD+1
      LTESCO=LRMAX+1
      LCCMAX=LTESCO+36
      LEL0=LCCMAX+1
      LH=LEL0+1
      LHMIN=LH+1
      LHMXI=LHMIN+1
      LHU=LHMXI+1
      LRC=LHU+1
      LTN=LRC+1
      LUROUND=LTN+1
      LALPHA=LUROUND+1
      LX2=LALPHA+1
      LT0=LX2+1
      LTLAST=LT0+1
      LTOUTC=LTLAST+1
      LCNTL=LTOUTC+1
      LRINFO=LCNTL+10

      LMBAND=3
      LNE=4
      LCHOUT=8
      LINDEX=19
      LVALUE=20
      LILLIN=21
      LINIT=LILLIN+1
      LLYH=LINIT+1
      LLEWT=LLYH+1
      LLACOR=LLEWT+1
      LLSAVR=LLACOR+1
      LLWM=LLSAVR+1
      LLIWM=LLWM+1
      LMXSTEP=LLIWM+1
      LMXHNIL=LMXSTEP+1
      LNHNIL=LMXHNIL+1
      LNTREP=LNHNIL+1
      LNSLAST=LNTREP+1
      LNYH=LNSLAST+1
      LIALTH=LNYH+1
      LIPUP=LIALTH+1
      LLMAX=LIPUP+1
      LMEO=LLMAX+1
      LNQNYH=LMEO+1
      LNSLP=LNQNYH+1
      LICF=LNSLP+1
      LIERPJ=LICF+1
      LIERSL=LIERPJ+1
      LJCUR=LIERSL+1
      LJSTART=LJCUR+1
      LKFLAG=LJSTART+1
      LL=LKFLAG+1
      LMETH=LL+1
      LMITER=LMETH+1
      LMAXORD=LMITER+1
      LMAXCOR=LMAXORD+1
      LMSBP=LMAXCOR+1
      LMXNCF=LMSBP+1
      LN=LMXNCF+1
      LNQ=LN+1
      LNST=LNQ+1
      LNRE=LNST+1
      LNJE=LNRE+1
      LNQU=LNJE+1
      LLG0=LNQU+1
      LLG1=LLG0+1
      LLGX=LLG1+1
      LIMAX=LLGX+1
      LLAST=LIMAX+1
      LIRFND=LLAST+1
      LITASKC=LIRFND+1
      LNGC=LITASKC+1
      LNGE=LNGC+1
      LMESFLG=LNGE+1
      LLUNIT=LMESFLG+1
      LKEEP=LLUNIT+1
      LINFO=LKEEP+20
      LICNTL=LINFO+40
      LIDX=LICNTL+20
C------------------------------------------------------------------------

C-----------------------------------------------------------------------
C BLOCK A.
C THIS CODE BLOCK IS EXECUTED ON EVERY CALL.
C IT TESTS ISTATE AND ITASK FOR LEGALITY AND BRANCHES APPROPIATELY.
C IF ISTATE .GT. 1 BUT THE FLAG INIT SHOWS THAT INITIALIZATION HAS
C NOT YET BEEN DONE, AN ERROR RETURN OCCURS.
C IF ISTATE = 0 OR 1 AND TOUT = T, JUMP TO BLOCK G AND RETURN
C IMMEDIATELY.
C-----------------------------------------------------------------------

      IF (IW(LMESFLG).LT.0.OR.IW(LMESFLG).GT.1) IW(LMESFLG)=0
      IF (IW(LLUNIT).GT.0) IW(LMESFLG)=1
      IF (ISTATE.LT.0.OR.ISTATE.GT.3) GOTO 601
      IF (ITASK.LT.1.OR.ITASK.GT.5) GOTO 602
      IW(LITASKC)=ITASK
      IF (ISTATE.LE.1) GOTO 10
      IF (IW(LINIT).EQ.0) GOTO 603
      IF (ISTATE.EQ.2) GOTO 200
      GOTO 20
 10   IW(LINIT)=0
      IF (TOUT.EQ.T) GO TO 430
 20   IW(LNTREP)=0

C-----------------------------------------------------------------------
C BLOCK B.
C THE NEXT CODE BLOCK IS EXECUTED FOR THE INITIAL CALL (ISTATE = 0 OR 1)
C OR FOR A CONTINUATION CALL WITH PARAMETER CHANGES (ISTATE = 3).
C IT CONTAINS CHECKING OF ALL INPUTS AND VARIOUS INITIALIZATIONS.
C
C FIRST CHECK LEGALITY OF THE NON-OPTIONAL INPUTS NEQ, ITOL, IOPT,
C MF, ML, AND MU.
C-----------------------------------------------------------------------

      IF (NEQ(1).LE.0) GOTO 604
      IF (ISTATE.LE.1) GOTO 25
      IF (NEQ(1).GT.IW(LN)) GOTO 605
 25   IW(LN)=NEQ(1)
      IF (ITOL.LT.1.OR.ITOL.GT.4) GOTO 606
      IF (IOPT.LT.0.OR.IOPT.GT.1) GOTO 607
      IW(LMETH)=MF/10
      IW(LMITER)=MF-10*IW(LMETH)
      IF (IW(LMETH).LT.1.OR.IW(LMETH).GT.2) GOTO 608
      IF (IW(LMITER).LE.0.OR.IW(LMITER).GT.8) GOTO 608
      IF (IW(LMITER).GE.5) GOTO 30
      IF (IW(LMITER).LT.3) GOTO 30
      ML=IW(1)
      MU=IW(2)
      IF (ML.LT.0.OR.ML.GE.IW(LN)) GOTO 609
      IF (MU.LT.0.OR.MU.GE.IW(LN)) GOTO 610
 30   CONTINUE
      IF (NG.LT.0) GOTO 630
      IF (ISTATE.LE.1) GOTO 35
      IF (IW(LIRFND).EQ.0.AND.NG.NE.IW(LNGC)) GOTO 631
 35   IW(LNGC)=NG

C NEXT PROCESS AND CHECK THE OPTIONAL INPUTS. --------------------------

      IW(LMSBP)=20
      IF (IOPT.EQ.1) GOTO 40
      IW(LMAXORD)=MORD(IW(LMETH))
      IW(LMXSTEP)=MXSTP0
      IW(LMXHNIL)=MXHNL0
      IF (ISTATE.LE.1) H0=0.0D0
      RW(LHMXI)=0.0D0
      RW(LHMIN)=0.0D0
      IW(9)=0
      GOTO 60
 40   IW(LMAXORD)=IW(5)
      IF (IW(LMAXORD).LT.0) GOTO 611
      IF (IW(LMAXORD).EQ.0) IW(LMAXORD)=100
      IW(LMAXORD)=MIN0(IW(LMAXORD),MORD(IW(LMETH)))
      IF (IW(LCHOUT).GT.0) IW(LMSBP)=IW(LCHOUT)
      IW(LMXSTEP)=IW(6)
      IF (IW(LMXSTEP).LT.0) GOTO 612
      IF (IW(LMXSTEP).EQ.0) IW(LMXSTEP)=MXSTP0
      IW(LMXHNIL)=IW(7)
      IF (IW(LMXHNIL).LT.0) GOTO 613
      IF (IW(LMXHNIL).EQ.0) IW(LMXHNIL)=MXHNL0
      IF (ISTATE.GT.1.AND.ISTATE.LT.4) GOTO 50
      H0=RW(5)
      IF ((TOUT-T)*H0.LT.0.0D0) GOTO 614
 50   HMAX=RW(6)
      IF (HMAX.LT.0.0D0) GOTO 615
      RW(LHMXI)=0.0D0
      IF (HMAX.GT.0.0D0) RW(LHMXI)=1.0D0/HMAX
      RW(LHMIN)=RW(7)
      IF (RW(LHMIN).LT.0.0D0) GOTO 616

C-----------------------------------------------------------------------
C SET WORK ARRAY POINTERS AND CHECK LENGTHS LRW AND LIW.
C POINTERS TO SEGMENTS OF RW AND IW ARE NAMED BY PREFIXING L TO
C THE NAME OF THE SEGMENT.  E.G., THE SEGMENT YH STARTS AT RW(LYH).
C SEGMENTS OF RW (IN ORDER) ARE DENOTED YH, WM, EWT, SAVR, ACOR.
C-----------------------------------------------------------------------

 60   CONTINUE
      IF (ISTATE.LE.1) IW(LNYH)=IW(LN)
      IW(LLG0)=21+254
      IW(LLG1)=IW(LLG0)+NG
      IW(LLGX)=IW(LLG1)+NG
      LYHNEW=IW(LLGX)+NG
      IF (ISTATE.LE.1) IW(LLYH)=LYHNEW
      IF (LYHNEW.EQ.IW(LLYH)) GOTO 62

C IF ISTATE = 3 AND NG WAS CHANGED, SHIFT YH TO ITS NEW LOCATION. ------

      LENYH=IW(LL)*IW(LNYH)
      IF (LRW.LT.LYHNEW-1+LENYH) GOTO 62
      I1=1
      IF (LYHNEW.GT.IW(LLYH)) I1=-1
      CALL DCOPY(LENYH,RW(IW(LLYH)),I1,RW(LYHNEW),I1)
      IW(LLYH)=LYHNEW
 62   CONTINUE
      IW(LLWM)=LYHNEW+(IW(LMAXORD)+1)*IW(LNYH)

CAS   Initialisation MA38

      CALL UMD2IN(IW(LICNTL),RW(LCNTL),IW(LKEEP))      
      IW(LICNTL+2)=0

CAS   Création index

      IF (IW(LMITER).LE.2) THEN
        CALL DENSIDX(IW(LIDX),NEQ,IW(LNE),IW(LINDEX),
     1               IW(LVALUE))
        IW(LMBAND)=NEQ(1)
      END IF
      IF (IW(LMITER).GE.3.AND.IW(LMITER).LE.4) THEN
        CALL BANDIDX(IW(LIDX),NEQ,IW(LNE),MU,ML,IW(LINDEX),
     1               IW(LVALUE))
        IW(LMBAND)=ML+MU+1
      END IF
      IF (IW(LMITER).GE.5.AND.IW(LMITER).LE.6) THEN
        CALL TBBIDX(IW(LIDX),NEQ,IW(LNE),IW(LINDEX),
     1              IW(LVALUE))
        IW(LMBAND)=1
      END IF
      IF (IW(LMITER).GE.7.AND.IW(LMITER).LE.8) THEN
        CALL SPAIDX(NEQ,IW(LNE),IW(LINDEX),
     1              IW(LVALUE))
        IW(LMBAND)=1
      END IF
      IF (IW(LMITER).GE.5.AND.IW(LMITER).LE.6) THEN
        IF (IW(1).LT.0) GOTO 609
      ENDIF

      LENWM=IW(LVALUE)+2+IW(LNE)+2*IW(LN)
      LENRW=IW(LLWM)+LENWM+3*IW(LN)-1
      IF(LENRW.GT.LRW) GOTO 617
      IW(LVALUE)=IW(LVALUE)+LRW-LENRW
      LENWM=IW(LVALUE)+2+IW(LNE)+2*IW(LN)
      IW(LLEWT)=IW(LLWM)+LENWM
      IW(LLSAVR)=IW(LLEWT)+IW(LN)
      IW(LLACOR)=IW(LLSAVR)+IW(LN)
      LENRW=LRW
      IW(17)=LENRW

      IW(LLIWM)=1
      LENIW=150+2*IW(LNE)+IW(LINDEX)
      IF(LENIW.GT.LIW) GOTO 618
      IW(LINDEX)=IW(LINDEX)+LIW-LENIW
      LENIW=LIW
      IW(18)=LENIW

C CHECK RTOL AND ATOL FOR LEGALITY. ------------------------------------

      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO I=1,IW(LN)
        IF (ITOL.GE.3) RTOLI=RTOL(I)
        IF (ITOL.EQ.2.OR.ITOL.EQ.4) ATOLI=ATOL(I)
        IF (RTOLI.LT.0.0D0) GOTO 619
        IF (ATOLI.LT.0.0D0) GOTO 620
      END DO
      IF (ISTATE.LE.1) GOTO 100 

C IF ISTATE = 3, SET FLAG TO SIGNAL PARAMETER CHANGES TO STODI. --------

      IW(LJSTART)=-1
      IF (IW(LNQ).LE.IW(LMAXORD)) GOTO 90

C MAXORD WAS REDUCED BELOW NQ.  COPY YH(*,MAXORD+2) INTO YDOTI.---------

      DO I=1,IW(LN)
        YDOTI(I)=RW(I+IW(LLWM)-1)
      END DO

C RELOAD WM(1) = RW(LWM), SINCE LWM MAY HAVE CHANGED. ---------------

 90   RW(IW(LLWM)) = DSQRT(RW(LUROUND))
      IF (IW(LN).EQ.IW(LNYH)) GOTO 200

C NEQ WAS REDUCED.  ZERO PART OF YH TO AVOID UNDEFINED REFERENCES. -----

      I1=IW(LLYH)+IW(LL)*IW(LNYH)
      I2=IW(LLYH)+(IW(LMAXORD)+1)*IW(LNYH)-1
      IF (I1.GT.I2) GOTO 200
      DO I=I1,I2
        RW(I)=0.0D0
      END DO
      GOTO 200

C-----------------------------------------------------------------------
C BLOCK C.
C THE NEXT BLOCK IS FOR THE INITIAL CALL ONLY (ISTATE = 0 OR 1).
C IT CONTAINS ALL REMAINING INITIALIZATIONS, THE CALL TO AINVG
C (IF ISTATE = 1), AND THE CALCULATION OF THE INITIAL STEP SIZE.
C THE ERROR WEIGHTS IN EWT ARE INVERTED AFTER BEING LOADED.
C-----------------------------------------------------------------------

 100  CONTINUE
      RW(LUROUND)=D1MACH()
      RW(LTN)=T
      IF (ITASK.NE.4.AND.ITASK.NE.5) GOTO 105
      TCRIT=RW(1)
      IF ((TCRIT-TOUT)*(TOUT-T).LT.0.0D0) GOTO 625
      IF (H0.NE.0.0D0.AND.(T+H0-TCRIT)*H0.GT.0.0D0) H0=TCRIT-T
 105  IW(LJSTART)=0
      RW(IW(LLWM))=DSQRT(RW(LUROUND))
      IW(LNHNIL)=0
      IW(LNST)=0
      IW(LNRE)=0
      IW(LNJE)=0
      IW(LNSLAST)=0
      RW(LHU)=0.0D0
      IW(LNQU)=0
      RW(LCCMAX)=0.3D0
      IF (IW(LMAXCOR).EQ.0) IW(LMAXCOR)=6
      IF (IW(LMXNCF).EQ.0) IW(LMXNCF)=5
      IW(LNQ)=1
      RW(LH)=1.0D0

C COMPUTE INITIAL DY/DT, IF NECESSARY, AND LOAD IT AND INITIAL Y INTO YH

      LYD0=IW(LLYH)+IW(LNYH)
      LP=IW(LLWM)+1

C CALCUL DES DERIVEES DES VARIABLES AU TEMPS INITIAL

      DO I=1,IW(LN)
        RW(I+IW(LLYH)-1)=Y(I)
        RW(I+LYD0-1)=YDOTI(I)
      END DO
      IF (ISTATE.NE.0) GOTO 133
      CALL EWSET (IW(LN),ITOL,RTOL,ATOL,RW(IW(LLYH)),
     1            RW(IW(LLEWT)))
      DO I=1,IW(LN)
        IF (RW(I+IW(LLEWT)-1).LE.0.0D0) GOTO 621
        RW(I+IW(LLEWT)-1)=1.0D0/RW(I+IW(LLEWT)-1)
      END DO
      HINI=H0
      IF (HINI.NE.0.0D0) GOTO 109
      TDIST=DABS(TOUT - T)
      W0=DMAX1(DABS(T), DABS(TOUT))
      IF (TDIST.LT.2.0D0*RW(LUROUND)*W0) GOTO 622
      HINI=0.001D0*TDIST
      SUM=VNORM(IW(LN),RW(LYD0),RW(IW(LLEWT)))
      IF (SUM.GT.0.5D0/HINI) HINI=0.5D0/SUM
      HINI=DSIGN(HINI,TOUT-T)
 109  RH=DABS(HINI)*RW(LHMXI)
      IF (RH.GT.1.0D0) HINI=HINI/RH
c     IW(LIRFND)=0
c     RW(LTOUTC)=TOUT
c     IF (IW(LNGC).NE.0) THEN
c       CALL SRCHEK(1,GEX,NEQ,IW(LN),Y,RW(IW(LLYH)),
c    1              IW(LNYH),RW(IW(LLG0)),RW(IW(LLG1)),
c    2              RW(IW(LLGX)),JROOT,IRT,YDOTI,RW(LTN),
c    3              RW(LH),RW(LHU),RW(LUROUND),RW(LT0),
c    4              RW(LTLAST),RW(LTOUTC),IW(LIRFND),
c    5              IW(LITASKC),IW(LL),IW(LNQ),IW(LNGC),
c    6              IW(LNGE),RPAR,LRP,IPAR,LIP,RW(LALPHA),
c    7              RW(LX2),IW(LIMAX),IW(LLAST),IW(LMESFLG),
c    8              IW(LLUNIT))
c       IF (IRT.NE.0) GOTO 632
c     END IF
c     RW(LH)=HINI

C DISCo MUST COMPUTE INITIAL DY/DT (LYD0 POINTS TO YH(*,2)). -----------

       CALL SAINVG(T,Y,YDOTI,NEQ,RES,ADDA,JAC,HINI,
     1             RW(IW(LLEWT)),IER,IW(LMITER),
     2             RW(IW(LLACOR)),RW(IW(LLYH)),
     3             IW(LNYH),RW(IW(LLSAVR)),
     4             RW(IW(LLWM)),IW(IW(LLIWM)),
     5             RW(LHMIN),RW(LUROUND),RPAR,LRP,IPAR,LIP,
     6             IW(LNRE),IW(LNJE),RW(LCNTL),RW(LRINFO),
     7             IW(LKEEP),IW(LINFO),IW(LICNTL))
      IW(LINIT)=1
      H0=HINI
      IF (IER) 560,110,565
 133  IF(ISTATE.EQ.1) GOTO 130
 110  CONTINUE
      RW(LTN)=T
      DO I=1,IW(LN)
        RW(I+IW(LLYH)-1)=Y(I)
        RW(I+LYD0-1)=YDOTI(I)
      END DO
c     IF (IW(LNGC).NE.0) THEN
c       CALL SRCHEK(3,GEX,NEQ,IW(LN),Y,RW(IW(LLYH)),
c    1              IW(LNYH),RW(IW(LLG0)),RW(IW(LLG1)),
c    2              RW(IW(LLGX)),JROOT,IRT,YDOTI,RW(LTN),
c    3              RW(LH),RW(LHU),RW(LUROUND),RW(LT0),
c    4              RW(LTLAST),RW(LTOUTC),IW(LIRFND),
c    5              IW(LITASKC),IW(LL),IW(LNQ),IW(LNGC),
c    6              IW(LNGE),RPAR,LRP,IPAR,LIP,RW(LALPHA),
c    7              RW(LX2),IW(LIMAX),IW(LLAST),IW(LMESFLG),
c    8              IW(LLUNIT))
c       IF (IRT.EQ.1) THEN
c         IW(LIRFND)=1
c         ISTATE=3
c         T=RW(LT0)
c         GOTO 425
c       END IF
c     END IF
 130  CONTINUE
      IW(LNQ)=1
      RW(LH)=1.0D0
      CALL EWSET(IW(LN),ITOL,RTOL,ATOL,RW(IW(LLYH)),
     1           RW(IW(LLEWT)))
      DO I=1,IW(LN)
        IF (RW(I+IW(LLEWT)-1).LE.0.0D0) GOTO 621
        RW(I+IW(LLEWT)-1)=1.0D0/RW(I+IW(LLEWT)-1)
      END DO

C-----------------------------------------------------------------------
C THE CODING BELOW COMPUTES THE STEP SIZE, H0, TO BE ATTEMPTED ON THE
C FIRST STEP, UNLESS THE USER HAS SUPPLIED A VALUE FOR THIS.
C FIRST CHECK THAT TOUT - T DIFFERS SIGNIFICANTLY FROM ZERO.
C A SCALAR TOLERANCE QUANTITY TOL IS COMPUTED, AS MAX(RTOL(I))
C IF THIS IS POSITIVE, OR MAX(ATOL(I)/ABS(Y(I))) OTHERWISE, ADJUSTED
C SO AS TO BE BETWEEN 100*UROUND AND 1.0E-3.
C THEN THE COMPUTED VALUE H0 IS GIVEN BY..
C                                      NEQ
C   H0**2 = TOL / ( W0**-2 + (1/NEQ) * SUM ( YDOT(I)/YWT(I) )**2  )
C                                       1
C WHERE   W0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(I) = I-TH COMPONENT OF INITIAL VALUE OF DY/DT,
C         YWT(I)  = EWT(I)/TOL  (A WEIGHT FOR Y(I)).
C THE SIGN OF H0 IS INFERRED FROM THE INITIAL VALUES OF TOUT AND T.
C-----------------------------------------------------------------------

      IF (H0.NE.0.0D0) GOTO 180
      TDIST=DABS(TOUT-T)
      W0=DMAX1(DABS(T),DABS(TOUT))
      IF (TDIST.LT.2.0D0*RW(LUROUND)*W0) GOTO 622
      TOL=RTOL(1)
      IF (ITOL.LE.2) GOTO 145
      DO I=1,IW(LN)
        TOL=DMAX1(TOL,RTOL(I))
      END DO
 145  IF (TOL.GT.0.0D0) GOTO 160
      ATOLI=ATOL(1)
      DO I=1,IW(LN)
        IF (ITOL.EQ.2.OR.ITOL.EQ.4) ATOLI=ATOL(I)
        AYI=DABS(Y(I))
        IF (AYI.NE.0.0D0) TOL=DMAX1(TOL,ATOLI/AYI)
      END DO
 160  TOL=DMAX1(TOL,100.0D0*RW(LUROUND))
      TOL=DMIN1(TOL,0.001D0)
      SUM=VNORM(IW(LN),RW(LYD0),RW(IW(LLEWT)))
      SUM=1.0D0/(TOL*W0*W0)+TOL*SUM**2
      H0=1.0D0/DSQRT(SUM)
      H0=DMIN1(H0,TDIST)
      H0=DSIGN(H0,TOUT-T)

C ADJUST H0 IF NECESSARY TO MEET HMAX BOUND. ---------------------------

 180  RH=DABS(H0)*RW(LHMXI)
      IF (RH.GT.1.0D0) H0=H0/RH

C LOAD H WITH H0 AND SCALE YH(*,2) BY H0. ------------------------------

      RW(LH)=H0
      DO I=1,IW(LN)
        RW(I+LYD0-1)=H0*RW(I+LYD0-1)
      END DO

C CHECK FOR A ZERO OF G AT T. ------------------------------------------

      IW(LIRFND)=0
      RW(LTOUTC)=TOUT
      IF (IW(LNGC).NE.0) THEN
        CALL SRCHEK(1,GEX,NEQ,IW(LN),Y,RW(IW(LLYH)),
     1              IW(LNYH),RW(IW(LLG0)),RW(IW(LLG1)),
     2              RW(IW(LLGX)),JROOT,IRT,YDOTI,RW(LTN),
     3              RW(LH),RW(LHU),RW(LUROUND),RW(LT0),
     4              RW(LTLAST),RW(LTOUTC),IW(LIRFND),
     5              IW(LITASKC),IW(LL),IW(LNQ),IW(LNGC),
     6              IW(LNGE),RPAR,LRP,IPAR,LIP,RW(LALPHA),
     7              RW(LX2),IW(LIMAX),IW(LLAST),IW(LMESFLG),
     8              IW(LLUNIT))
        IF (IRT.NE.0) GOTO 632
      END IF
      GOTO 270
C-----------------------------------------------------------------------
C BLOCK D.
C THE NEXT CODE BLOCK IS FOR CONTINUATION CALLS ONLY (ISTATE = 2 OR 3)
C AND IS TO CHECK STOP CONDITIONS BEFORE TAKING A STEP.
C FIRST, RCHEK IS CALLED TO CHECK FOR A ROOT WITHIN THE LAST STEP
C TAKEN, OTHER THAN THE LAST ROOT FOUND THERE, IF ANY.
C IF ITASK = 2 OR 5, AND Y(TN) HAS NOT YET BEEN RETURNED TO THE USER
C BECAUSE OF AN INTERVENING ROOT, RETURN THROUGH BLOCK G.
C-----------------------------------------------------------------------

 200  IW(LNSLAST)=IW(LNST)
      IRFP=IW(LIRFND)
      IF (IW(LNGC).NE.0) THEN
        IF (ITASK.EQ.1.OR.ITASK.EQ.4) RW(LTOUTC)=TOUT
        CALL SRCHEK(2,GEX,NEQ,IW(LN),Y,RW(IW(LLYH)),
     1              IW(LNYH),RW(IW(LLG0)),RW(IW(LLG1)),
     2              RW(IW(LLGX)),JROOT,IRT,YDOTI,RW(LTN),
     3              RW(LH),RW(LHU),RW(LUROUND),RW(LT0),
     4              RW(LTLAST),RW(LTOUTC),IW(LIRFND),
     5              IW(LITASKC),IW(LL),IW(LNQ),IW(LNGC),
     6              IW(LNGE),RPAR,LRP,IPAR,LIP,RW(LALPHA),
     7              RW(LX2),IW(LIMAX),IW(LLAST),IW(LMESFLG),
     8              IW(LLUNIT))
        IF (IRT.EQ.1) THEN
          IW(LIRFND)=1
          ISTATE=3
          T=RW(LT0)
          GOTO 425
        END IF
      END IF
      IW(LIRFND)=0
      IF (IRFP.EQ.1.AND.RW(LTLAST).NE.RW(LTN).AND.ITASK.EQ.2) 
     1  GOTO 400
      GOTO (210,250,220,230,240), ITASK
 210  IF ((RW(LTN)-TOUT)*RW(LH).LT.0.0D0) GOTO 250
      CALL INTDY(IW(LN),TOUT,RW(LTN),0,RW(IW(LLYH)),
     1           IW(LNYH),Y,RW(LH),RW(LHU),RW(LUROUND),
     2           IW(LL),IW(LNQ),IFLAG,IW(LMESFLG),
     3           IW(LLUNIT))
      IF (IFLAG.NE.0) GOTO 627
      T=TOUT
      GOTO 420
 220  TP=RW(LTN)-RW(LHU)*(1.0D0+100.0D0*RW(LUROUND))
      IF ((TP-TOUT)*RW(LH).GT.0.0D0) GOTO 623
      IF ((RW(LTN)-TOUT)*RW(LH).LT.0.0D0) GOTO 250
      T=RW(LTN)
      GOTO 400
 230  TCRIT=RW(1)
      IF ((RW(LTN)-TCRIT)*RW(LH).GT.0.0D0) GOTO 624
      IF ((TCRIT-TOUT)*RW(LH).LT.0.0D0) GOTO 625
      IF ((RW(LTN)-TOUT)*RW(LH).LT.0.0D0) GOTO 245
      CALL INTDY(IW(LN),TOUT,RW(LTN),0,RW(IW(LLYH)),
     1           IW(LNYH),Y,RW(LH),RW(LHU),RW(LUROUND),
     2           IW(LL),IW(LNQ),IFLAG,IW(LMESFLG),
     3           IW(LLUNIT))
      IF (IFLAG.NE.0) GOTO 627
      T=TOUT
      GOTO 420
 240  TCRIT=RW(1)
      IF ((RW(LTN)-TCRIT)*RW(LH).GT.0.0D0) GOTO 624
 245  HMX=DABS(RW(LTN))+DABS(RW(LH))
      IHIT=DABS(RW(LTN)-TCRIT).LE.100.0D0*RW(LUROUND)*HMX
      IF(IHIT) T=TCRIT
      IF (IRFP.EQ.1.AND.RW(LTLAST).NE.RW(LTN).AND.ITASK.EQ.5) 
     1  GOTO 400
      IF (IHIT) GOTO 400
      TNEXT=RW(LTN)+RW(LH)*(1.0D0+4.0D0*RW(LUROUND))
      IF ((TNEXT-TCRIT)*RW(LH).LE.0.0D0) GOTO 250
      RW(LH)=(TCRIT-RW(LTN))*(1.0D0-4.0D0*RW(LUROUND))

CAS   CORRECTION DU BUG INDIQUE PAR JLD 26/06/97

      IF (ISTATE.EQ.2.AND.IW(LJSTART).GE.0) IW(LJSTART)=-2

C-----------------------------------------------------------------------
C BLOCK E.
C THE NEXT BLOCK IS NORMALLY EXECUTED FOR ALL CALLS AND CONTAINS
C THE CALL TO THE ONE-STEP CORE INTEGRATOR STODI.
C
C THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
C
C FIRST CHECK FOR TOO MANY STEPS BEING TAKEN, UPDATE EWT (IF NOT AT
C START OF PROBLEM), CHECK FOR TOO MUCH ACCURACY BEING REQUESTED, AND
C CHECK FOR H BELOW THE ROUNDOFF LEVEL IN T.
C-----------------------------------------------------------------------

 250  CONTINUE
      IF ((IW(LNST)-IW(LNSLAST)).GE.IW(LMXSTEP)) GOTO 500
      CALL EWSET(IW(LN),ITOL,RTOL,ATOL,RW(IW(LLYH)),
     1           RW(IW(LLEWT)))
      DO I=1,IW(LN)
        IF (RW(I+IW(LLEWT)-1).LE.0.0D0) GOTO 510
        RW(I+IW(LLEWT)-1)=1.0D0/RW(I+IW(LLEWT)-1)
      END DO
 270  TOLSF=RW(LUROUND)*VNORM(IW(LN),RW(IW(LLYH)),
     1                           RW(IW(LLEWT)))
      IF (TOLSF.LE.1.0D0) GOTO 280
      TOLSF=TOLSF*2.0D0
      IF (IW(LNST).EQ.0) GOTO 626
      GOTO 520
 280  IF ((RW(LTN)+RW(LH)).NE.RW(LTN)) GOTO 290
      IW(LNHNIL)=IW(LNHNIL)+1
      IF (IW(LNHNIL).GT.IW(LMXHNIL)) GOTO 290
      CALL XERRWV("WARNING : Internal T (=R1) and H (=R2) are",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1 "          such that in the machine, T + H = T on the next step",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      CALL XERRWV(
     1 "          (H = step size). Solver will continue anyway",
     1   1,0,0,0,2,RW(LTN),RW(LH),IW(LMESFLG),IW(LLUNIT),0)
      IF (IW(LNHNIL).LT.IW(LMXHNIL)) GOTO 290
      CALL XERRWV("          Above warning has been issued I1 times.",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      CALL XERRWV(
     1 "          It will not be issued again for this problem",
     1   1,1,IW(LMXHNIL),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
 290  CONTINUE

C NOTE... SAVF IN STODI OCCUPIES THE SAME SPACE AS YDOTI IN DISCo.

      CALL SSTODI(NEQ,Y,RW(IW(LLYH)),IW(LNYH),
     1            RW(IW(LLYH)),RW(IW(LLEWT)),YDOTI,
     2            RW(IW(LLSAVR)),RW(IW(LLACOR)),
     3            RW(IW(LLWM)),IW(IW(LLIWM)),
     4            RW,IW,RW(LEL),RW(LELCO),RW(LTESCO),
     5            RES,ADDA,JAC,RPAR,LRP,IPAR,LIP,ALPHA)
      KGO=1-IW(LKFLAG)
      GOTO (300,530,540,400,550,560), KGO

C KGO = 1,SUCCESS. 2,ERROR TEST FAILURE. 3,CONVERGENCE FAILURE.
C       4,RES ORDERED RETURN. 5,RES RETURNED ERROR.

C-----------------------------------------------------------------------
C BLOCK F.
C THE FOLLOWING BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN FROM THE
C CORE INTEGRATOR (KFLAG = 0).  TEST FOR STOP CONDITIONS.
C THEN CALL RCHEK TO CHECK FOR A ROOT WITHIN THE LAST STEP.
C THEN, IF NO ROOT WAS FOUND, CHECK FOR STOP CONDITIONS.
C-----------------------------------------------------------------------

 300  IW(LINIT)=1
      IF (IW(LNGC).NE.0) THEN
        CALL SRCHEK(3,GEX,NEQ,IW(LN),Y,RW(IW(LLYH)),
     1              IW(LNYH),RW(IW(LLG0)),RW(IW(LLG1)),
     2              RW(IW(LLGX)),JROOT,IRT,YDOTI,RW(LTN),
     3              RW(LH),RW(LHU),RW(LUROUND),RW(LT0),
     4              RW(LTLAST),RW(LTOUTC),IW(LIRFND),
     5              IW(LITASKC),IW(LL),IW(LNQ),IW(LNGC),
     6              IW(LNGE),RPAR,LRP,IPAR,LIP,RW(LALPHA),
     7              RW(LX2),IW(LIMAX),IW(LLAST),IW(LMESFLG),
     8              IW(LLUNIT))
        IF (IRT.EQ.1) THEN
          IW(LIRFND)=1
          ISTATE=3
          T=RW(LT0)
          GOTO 425
        END IF
      END IF
      GOTO (310,400,330,340,350), ITASK

C ITASK = 1.  IF TOUT HAS BEEN REACHED, INTERPOLATE. -------------------

 310  IF ((RW(LTN)-TOUT)*RW(LH).LT.0.0D0) GOTO 250
      CALL INTDY(IW(LN),TOUT,RW(LTN),0,RW(IW(LLYH)),
     1           IW(LNYH),Y,RW(LH),RW(LHU),RW(LUROUND),
     2           IW(LL),IW(LNQ),IFLAG,IW(LMESFLG),
     3           IW(LLUNIT))
      T=TOUT
      GOTO 420

C ITASK = 3.  JUMP TO EXIT IF TOUT WAS REACHED. ------------------------

 330  IF ((RW(LTN)-TOUT)*RW(LH).GE.0.0D0) GOTO 400
      GOTO 250

C ITASK = 4.  SEE IF TOUT OR TCRIT WAS REACHED.  ADJUST H IF NECESSARY.

 340  IF ((RW(LTN)-TOUT)*RW(LH).LT.0.0D0) GOTO 345
      CALL INTDY(IW(LN),TOUT,RW(LTN),0,RW(IW(LLYH)),
     1           IW(LNYH),Y,RW(LH),RW(LHU),RW(LUROUND),
     2           IW(LL),IW(LNQ),IFLAG,IW(LMESFLG),
     3           IW(LLUNIT))
      T=TOUT
      GOTO 420
 345  HMX=DABS(RW(LTN))+DABS(RW(LH))
      IHIT=DABS(RW(LTN)-TCRIT).LE.100.0D0*RW(LUROUND)*HMX
      IF (IHIT) GOTO 400
      TNEXT=RW(LTN)+RW(LH)*(1.0D0+4.0D0*RW(LUROUND))
      IF ((TNEXT-TCRIT)*RW(LH).LE.0.0D0) GOTO 250
      RW(LH)=(TCRIT-RW(LTN))*(1.0D0-4.0D0*RW(LUROUND))

CAS   CORRECTION DU BUG INDIQUE PAR JLD 26/06/97

      IF (IW(LJSTART).GE.0) IW(LJSTART)=-2
      GOTO 250

C ITASK = 5.  SEE IF TCRIT WAS REACHED AND JUMP TO EXIT. ---------------

 350  HMX=DABS(RW(LTN))+DABS(RW(LH))
      IHIT=DABS(RW(LTN)-TCRIT).LE.100.0D0*RW(LUROUND)*HMX

C-----------------------------------------------------------------------
C BLOCK G.
C THE FOLLOWING BLOCK HANDLES ALL SUCCESSFUL RETURNS FROM DISCo.
C IF ITASK .NE. 1, Y IS LOADED FROM YH AND T IS SET ACCORDINGLY.
C ISTATE IS SET TO 2, THE ILLEGAL INPUT COUNTER IS ZEROED, AND THE
C OPTIONAL OUTPUTS ARE LOADED INTO THE WORK ARRAYS BEFORE RETURNING.  IF
C ISTATE = 0 OR 1 AND TOUT = T, THERE IS A RETURN WITH NO ACTION TAKEN,
C EXCEPT THAT IF THIS HAS HAPPENED REPEATEDLY, THE RUN IS TERMINATED.
C-----------------------------------------------------------------------

 400  DO I=1,IW(LN)
        Y(I)=RW(I+IW(LLYH)-1)
      END DO
      T=RW(LTN)
      IF (ITASK.NE.4.AND.ITASK.NE.5) GOTO 420
      IF (IHIT) T=TCRIT
 420  ISTATE=2
 425  CONTINUE
      IF (IW(LKFLAG).EQ.-3) ISTATE=3
      IW(LILLIN)=0
      RW(11)=RW(LHU)
      RW(12)=RW(LH)
      RW(13)=RW(LTN)
      IW(11)=IW(LNST)
      IW(12)=IW(LNRE)
      IW(13)=IW(LNJE)
      IW(14)=IW(LNQU)
      IW(15)=IW(LNQ)
      IW(10)=IW(LNGE)
      RW(LTLAST)=T
      RETURN

 430  IW(LNTREP)=IW(LNTREP)+1
      IF (IW(LNTREP).LT.5) RETURN
      CALL XERRWV(
     1  "ERROR :   Repeated calls with ISTATE= 0 or 1 and TOUT= T(=R1)",
     1   1,0,0,0,1,T,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      GO TO 800

C-----------------------------------------------------------------------
C BLOCK H.
C THE FOLLOWING BLOCK HANDLES ALL UNSUCCESSFUL RETURNS OTHER THAN
C THOSE FOR ILLEGAL INPUT.  FIRST THE ERROR MESSAGE ROUTINE IS CALLED.
C IF THERE WAS AN ERROR TEST OR CONVERGENCE TEST FAILURE, IMXER IS SET.
C THEN Y IS LOADED FROM YH, T IS SET TO TN, AND THE ILLEGAL INPUT
C COUNTER ILLIN IS SET TO 0.  THE OPTIONAL OUTPUTS ARE LOADED INTO
C THE WORK ARRAYS BEFORE RETURNING.
C-----------------------------------------------------------------------

C THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE REACHING TOUT. ----------

 500  CALL XERRWV("ERROR :   At current T (=R1), MXSTEP (=I1) steps",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          taken on this call before reaching TOUT",
     1   1,1,IW(LMXSTEP),0,1,RW(LTN),0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-100
      GOTO 580

C EWT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM). ----------------

 510  EWTI=RW(IW(LLEWT)+I-1)
      CALL XERRWV("ERROR :   At T (=R1), EWT(I1) has become R2 <= 0.",
     1   1,1,I,0,2,RW(LTN),EWTI,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-600
      GOTO 590

C TOO MUCH ACCURACY REQUESTED FOR MACHINE PRECISION. -------------------

 520  CALL XERRWV("ERROR :   At T (=R1), too much accuracy requested",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          for precision of machine. See TOLSF (=R2)",
     1   1,0,0,0,2,RW(LTN),TOLSF,IW(LMESFLG),IW(LLUNIT),0)
      RW(14)=TOLSF
      ISTATE=-200
      GOTO 590

C KFLAG = -1.  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H) = HMIN. -----

 530  CALL XERRWV("ERROR :   At T (=R1) and step size H (=R2), the",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1 "          error test failed repeatedly or with ABS(H) = HMIN",
     1   1,0,0,0,2,RW(LTN),RW(LH),IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-400
      GOTO 570

C KFLAG = -2.  CONVERGENCE FAILED REPEATEDLY OR WITH ABS(H) = HMIN. ----

 540  CALL XERRWV("ERROR :   At T (=R1) and step size H (=R2), the",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          corrector convergence failed repeatedly",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      CALL XERRWV("          or with ABS(H) = HMIN",
     1   1,0,0,0,2,RW(LTN),RW(LH),IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-500
      GOTO 570

C IRES = 3 RETURNED BY RES, DESPITE RETRIES BY STODI. ------------------

 550  CALL XERRWV("ERROR :   At T (=R1) residual routine returned",
     1   1,0,0,0,0,0.0D0,0.0D0 ,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          error IRES = 3 repeatedly",1,
     1   0,0,0,1,RW(LTN),0.0D0 ,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-700
      GOTO 590

 560  IF (IW(LINFO).EQ.-3) THEN
        LENIW=LENIW-IW(LINDEX)+IW(LINFO+18)
        IW(18)=LENIW
        GOTO 618
      END IF 
      IF ((IW(LINFO).EQ.-4).OR.(IW(LINFO).EQ.-5)) THEN
        LENRW=LENRW-IW(LVALUE)+IW(LINFO+20)
        IW(17)=LENRW
        GOTO 617
      END IF 
      CALL XERRWV(
     1  "ERROR :   SAINVG or SSTODI failed to perform one step",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1  "          UMD2FA, UMD2RF or UMD2SO returned INFO(1) = (I1)",
     2   1,1,IW(LINFO),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-800
      RETURN


 565  CALL XERRWV("ERROR :   attempt to initialize dX/dt failed",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1 "          because residual routine set its error flag",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      CALL XERRWV("          to IRES = (I1)",
     1   1,1,IER,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-801
      RETURN

C COMPUTE IMXER IF RELEVANT. -------------------------------------------

 570  BIG=0.0D0
      IMXER=1
      DO I=1,IW(LN)
        SIZE=DABS(RW(I+IW(LLACOR)-1)*RW(I+IW(LLEWT)-1))
        IF (BIG.LT.SIZE) THEN
          BIG=SIZE
          IMXER=I
        END IF
      END DO
      IW(16)=IMXER

C COMPUTE RESIDUAL IF RELEVANT. ----------------------------------------

 580  LYD0=IW(LLYH)+IW(LNYH)
      DO I=1,IW(LN)
        RW(I+IW(LLSAVR)-1)=RW(I+LYD0-1)/RW(LH)
        Y(I)=RW(I+IW(LLYH)-1)
      END DO
      IRES=1
      IF(IW(LMITER).GE.5) THEN
        INDPRP=0
        CON=0.
        CALL RES(NEQ,RW(LTN),Y,RW(IW(LLSAVR)),YDOTI,
     1           RW(IW(LLWM)+2),RPAR,LRP,IPAR,LIP,INDPRP,CON,
     2           IER,IW(IW(LLIWM)+20))
        ELSE
        CALL RES(NEQ,RW(LTN),Y,RW(IW(LLSAVR)),YDOTI,RPAR,
     1           LRP,IPAR,LIP)
      END IF
      IW(LNRE)=IW(LNRE)+1
      IF (IRES.LE.1) GOTO 595
      CALL XERRWV("ERROR :   Residual routine set its flag IRES",
     1    1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          to (I1) when called for final output.",
     1     1,1,IRES,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      GOTO 595

C SET Y VECTOR, T, ILLIN, AND OPTIONAL OUTPUTS. ------------------------

 590  DO I=1,IW(LN)
        Y(I)=RW(I+IW(LLYH)-1)
      END DO
 595  T=RW(LTN)
      IW(LILLIN)=0
      RW(11)=RW(LHU)
      RW(12)=RW(LH)
      RW(13)=RW(LTN)
      IW(11)=IW(LNST)
      IW(12)=IW(LNRE)
      IW(13)=IW(LNJE)
      IW(14)=IW(LNQU)
      IW(15)=IW(LNQ)
      IW(10)=IW(LNGE)
      RW(LTLAST)=T
      RETURN

C-----------------------------------------------------------------------
C BLOCK I.
C THE FOLLOWING BLOCK HANDLES ALL ERROR RETURNS DUE TO ILLEGAL INPUT
C (ISTATE = -30x), AS DETECTED BEFORE CALLING THE CORE INTEGRATOR.
C FIRST THE ERROR MESSAGE ROUTINE IS CALLED.  THEN IF THERE HAVE BEEN
C 5 CONSECUTIVE SUCH RETURNS JUST BEFORE THIS CALL TO THE SOLVER,
C THE RUN IS HALTED.
C-----------------------------------------------------------------------

 601  CALL XERRWV("ERROR :   ISTATE (=I1) illegal",
     1   1,1,ISTATE,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-300
      GOTO 700
 602  CALL XERRWV("ERROR :   ITASK (=I1) illegal",
     1   1,1,ITASK,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-301
      GOTO 700
 603  CALL XERRWV("ERROR :   ISTATE > 1 but DISCo not initialized",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-302
      GOTO 700
 604  CALL XERRWV("ERROR :   NEQ (=I1) < 1 illegal",
     1   1,1,NEQ(1),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-303
      GOTO 700
 605  CALL XERRWV("ERROR :   ISTATE = 3 and NEQ increased (I1 to I2)",
     1   1,2,IW(LN),NEQ(1),0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-304
      GOTO 700
 606  CALL XERRWV("ERROR :   ITOL (=I1) illegal",
     1   1,1,ITOL,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-305
      GOTO 700
 607  CALL XERRWV("ERROR :   IOPT (=I1) illegal",
     1   1,1,IOPT,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-306
      GOTO 700
 608  CALL XERRWV("ERROR :   MF (=I1) illegal",
     1   1,1,MF,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-307
      GOTO 700
 609  CALL XERRWV("ERROR :   ML(=I1) illegal, < 0 or >= NEQ(=I2)",
     1   1,2,ML,NEQ(1),0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-308
      GOTO 700
 610  CALL XERRWV("ERROR :   MU(=I1) illegal, < 0 or >= NEQ(=I2)",
     1   1,2,MU,NEQ(1),0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-309
      GOTO 700
 611  CALL XERRWV("ERROR :   MAXORD (=I1) < 0 illegal",
     1   1,1,IW(LMAXORD),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-310
      GOTO 700
 612  CALL XERRWV("ERROR :   MXSTEP (=I1) < 0 illegal",
     1   1,1,IW(LMXSTEP),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-311
      GOTO 700
 613  CALL XERRWV("ERROR :   MXHNIL (=I1) < 0 illegal",
     1   1,1,IW(LMXHNIL),0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-312
      GOTO 700
 614  CALL XERRWV("ERROR :  TOUT (=R1) behind T (=R2)",
     1   1,0,0,0,2,TOUT,T,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("         integration direction is given by H0 (=R1)",
     1   1,0,0,0,1,H0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-313
      GOTO 700
 615  CALL XERRWV("ERROR :   HMAX (=R1) < 0.0 illegal",
     1   1,0,0,0,1,HMAX,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-314
      GOTO 700
 616  CALL XERRWV("ERROR :   HMIN (=R1) < 0.0 illegal",
     1   1,0,0,0,1,RW(LHMIN),0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-315
      GOTO 700
 617  CALL XERRWV(
     1  "ERROR :   RWORK lenght needed, LENRW (=I1), exceeds LRW (=I2)",
     1   1,2,LENRW,LRW,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-316
      GOTO 700
 618  CALL XERRWV(
     1  "ERROR :   IWORK lenght needed, LENIW (=I1), exceeds LIW (=I2)",
     1   1,2,LENIW,LIW,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-317
      GOTO 700
 619  CALL XERRWV("ERROR :   RTOL(=I1) is R1 < 0.0 illegal",
     1   1,1,I,0,1,RTOLI,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-318
      GOTO 700
 620  CALL XERRWV("ERROR :   ATOL(=I1) is R1 < 0.0 illegal",
     1   1,1,I,0,1,ATOLI,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-319
      GOTO 700
 621  EWTI=RW(IW(LLEWT)+I-1)
      CALL XERRWV("ERROR :   EWT(=I1) is R1 < 0.0 illegal",
     1   1,1,I,0,1,EWTI, 0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-320
      GOTO 700
 622  CALL XERRWV(
     1  "ERROR :   TOUT (=R1) too close to T(=R2) to start integration",
     1   1,0,0,0,2,TOUT,T,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-321
      GOTO 700
 623  CALL XERRWV(
     1  "ERROR :   ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)",
     1   1,1,ITASK,0,2,TOUT,TP,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-322
      GOTO 700
 624  CALL XERRWV(
     1  "ERROR :   ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)",
     1   1,0,0,0,2,TCRIT,RW(LTN),IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-323
      GOTO 700
 625  CALL XERRWV(
     1  "ERROR :   ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)",
     1   1,0,0,0,2,TCRIT,TOUT,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-324
      GOTO 700
 626  CALL XERRWV("ERROR :   At start of problem, too much accuracy",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1  "          requested for precision of machine. See TOLSF (=R1)",
     1   1,0,0,0,1,TOLSF,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      RW(14)=TOLSF
      ISTATE=-325
      GOTO 700
 627  CALL XERRWV("ERROR :   Trouble from INTDY. ITASK = I1, TOUT = R1",
     1   1,1,ITASK,0,1,TOUT,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-326
      GOTO 700
 630  CALL XERRWV("ERROR :   NG (=I1) < 0 illegal",
     1   1,1,NG,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      ISTATE=-327
      GOTO 700
 631  CALL XERRWV("ERROR :   NG changed (from I1 to I2) illegally,",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV(
     1 "          i.e. not immediately after a root was found",
     1   1,2,IW(LNGC),NG,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-328
      GOTO 700
 632  CALL XERRWV("ERROR :   One or more components of G has a root",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
      CALL XERRWV("          too near to the initial point",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      ISTATE=-329

 700  IF (IW(LILLIN).EQ.5) GOTO 710
      IW(LILLIN)=IW(LILLIN)+1
      RETURN
 710  CALL XERRWV("ERROR :   Repeated occurrences of illegal input",
     1   1,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),1)
 800  CALL XERRWV("          run aborted. Apparent infinite loop",
     1   2,0,0,0,0,0.0D0,0.0D0,IW(LMESFLG),IW(LLUNIT),0)
      RETURN
C----------------------- Fin de la routine DISCo -----------------------
      END
      SUBROUTINE EWSET(N,ITOL,RTOL,ATOL,YCUR,EWT)
                         
C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 15 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------
C THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR EWT ACCORDING TO                 
C     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I),  I = 1,...,N,                    
C WITH THE SUBSCRIPT ON RTOL AND/OR ATOL POSSIBLY REPLACED BY 1 ABOVE,          
C DEPENDING ON THE VALUE OF ITOL.                                               
C-----------------------------------------------------------------------        
CLLL. OPTIMIZE                                                                  

      IMPLICIT NONE

      INTEGER N,ITOL,I
      DOUBLE PRECISION RTOL(*),ATOL(*),YCUR(*),EWT(*),ATOLI,RTOLI

      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO I=1,N
        IF (ITOL.GE.3) RTOLI=RTOL(I)
        IF (ITOL.EQ.2.OR.ITOL.EQ.4) ATOLI=ATOL(I)
        EWT(I)=RTOLI*DABS(YCUR(I))+ATOLI
      END DO
      RETURN
C----------------------- Fin de la routine EWSET -----------------------        
      END
C---------------------------------------------------------------------
C     Routines de Manipulation de Vecteurs
C---------------------------------------------------------------------
      SUBROUTINE ICOPY(N,SX,INCX,SY,INCY)
C                                                                               
C     COPIES A VECTOR, X, TO A VECTOR, Y.                                       
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.                            
C     JACK DONGARRA, LINPACK, 3/11/78.                                          
C                                                                               
      IMPLICIT NONE
      INTEGER SX(*),SY(*)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C                                                                               
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C                                                                               
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                        
C          NOT EQUAL TO 1                                                       
C                                                                               
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C                                                                               
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                                    
C                                                                               
C                                                                               
C        CLEAN-UP LOOP                                                          
C                                                                               
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE INTDY(N,T,TN,K,YH,NYH,DKY,H,HU,UROUND,L,NQ,
     1                 IFLAG,MESFLG,LUNIT)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      INTEGER K,NYH,IFLAG,L,NQ,I,IC,J,JB,JB2,JJ,JJ1,JP1,N,
     1        MESFLG,LUNIT
      DOUBLE PRECISION T,YH(NYH,*),DKY(*),H,HU,TN,UROUND,C,R,S,TP

C-----------------------------------------------------------------------        
C INTDY COMPUTES INTERPOLATED VALUES OF THE K-TH DERIVATIVE OF THE              
C DEPENDENT VARIABLE VECTOR Y, AND STORES IT IN DKY.  THIS ROUTINE              
C IS CALLED WITHIN THE PACKAGE WITH K = 0 AND T = TOUT, BUT MAY                 
C ALSO BE CALLED BY THE USER FOR ANY K UP TO THE CURRENT ORDER.                 
C (SEE DETAILED INSTRUCTIONS IN THE USAGE DOCUMENTATION.)                       
C-----------------------------------------------------------------------        
C THE COMPUTED VALUES IN DKY ARE GOTTEN BY INTERPOLATION USING THE              
C NORDSIECK HISTORY ARRAY YH.  THIS ARRAY CORRESPONDS UNIQUELY TO A             
C VECTOR-VALUED POLYNOMIAL OF DEGREE NQCUR OR LESS, AND DKY IS SET              
C TO THE K-TH DERIVATIVE OF THIS POLYNOMIAL AT T.                               
C THE FORMULA FOR DKY IS..                                                      
C              Q                                                                
C  DKY(I)  =  SUM  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)               
C             J=K                                                               
C WHERE  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.          
C THE QUANTITIES  NQ = NQCUR, L = NQ+1, N = NEQ, TN, AND H ARE                  
C COMMUNICATED BY COMMON.  THE ABOVE SUM IS DONE IN REVERSE ORDER.              
C IFLAG IS RETURNED NEGATIVE IF EITHER K OR T IS OUT OF BOUNDS.                 
C-----------------------------------------------------------------------        

      IFLAG=0
      IF (K.LT.0.OR.K.GT.NQ) GOTO 80
      TP=TN-HU*(1.0D0+100.0D0*UROUND)
      IF ((T-TP)*(T-TN).GT.0.0D0) GOTO 90
      S=(T-TN)/H
      IC=1
      IF (K.NE.0) THEN
        JJ1=L-K
        DO JJ=JJ1,NQ
          IC=IC*JJ
        END DO
      END IF
      C=DFLOAT(IC)
      DO I=1,N
        DKY(I)=C*YH(I,L)
      END DO
      IF (K.NE.NQ) THEN
        JB2=NQ-K
        DO JB=1,JB2
          J=NQ-JB
          JP1=J+1
          IC=1
          IF (K.NE.0) THEN
            JJ1=JP1-K
            DO JJ=JJ1,J
              IC=IC*JJ
            END DO
          END IF
          C=DFLOAT(IC)
          DO I=1,N
            DKY(I)=C*YH(I,JP1)+S*DKY(I)
          END DO
        END DO
        IF (K.EQ.0) RETURN
      END IF
      R=H**(-K)
      DO I=1,N
        DKY(I)=R*DKY(I)
      END DO
      RETURN
 80   CALL XERRWV("ERROR :  (INTDY)  K (=I1) illegal",
     1   1,1,K,0,0,0.0D0,0.0D0,MESFLG,LUNIT,1)
      IFLAG=-1
      RETURN
 90   CALL XERRWV("ERROR :  (INTDY)  T (=R1) illegal",
     1   1,0,0,0,1,T,0.0D0,MESFLG,LUNIT,1)
      CALL XERRWV(
     1  "          T not in interval TCUR - HU (= R1) to TCUR (=R2)",
     1   1,0,0,0,2,TP,TN,MESFLG,LUNIT,0)
      IFLAG=-2
      RETURN
C----------------------- Fin de la routine INTDY -----------------------        
      END
C-------------------------------------------------------------------
       SUBROUTINE NULLIDX(INDEX,LINDEX)
       
        IMPLICIT NONE
        INTEGER INDEX(*),I,LINDEX

        DO I=1,LINDEX
          INDEX(I)=0
        END DO
        RETURN
        END
      SUBROUTINE OPNUM(NEQ,Y,T,CON,EWT,RTEM,SAVR,S,WM,IWM,RES,RPAR,LRP,
     1                 IPAR,LIP,NRE,MITER)
C
      IMPLICIT NONE
      EXTERNAL RES
      INTEGER NEQ(*),IWM(*),MITER,IER,I,J,K,N,LENP,IRES,NRE,IPAR(*),
     1        LRP,LIP
      DOUBLE PRECISION Y(*),EWT(*),RTEM(*),SAVR(*),S(*),WM(*),YJ,SJ,T,
     1        CON,SRUR,R,FAC,RPAR(*)
C-----------------------------------------------------------------------
      LENP=IWM(4)
      DO I=1,LENP
        WM(2+I)=0.D0
      END DO
      N=NEQ(1)
      IRES=-1
      IF (MITER.NE.6) THEN
        CALL RES(NEQ,T,Y,S,SAVR,RPAR,LRP,IPAR,LIP)
      ELSE
        CALL RES(NEQ,T,Y,S,SAVR,WM(3),RPAR,LRP,IPAR,LIP,0,CON,IER,
     1           IWM(151))
      END IF
      NRE=NRE+1
      IF (IRES.GT.1) GOTO 600
      SRUR=WM(1)
      DO J=1,N
        YJ=Y(J)
        R=DMAX1(SRUR*DABS(YJ),0.01D0/EWT(J))
        Y(J)=Y(J)+R
        FAC=CON/R
        IF (MITER.NE.6) THEN
          CALL RES(NEQ,T,Y,S,RTEM,RPAR,LRP,IPAR,LIP)
        ELSE
          CALL RES(NEQ,T,Y,S,RTEM,WM(3),RPAR,LRP,IPAR,LIP,0,CON,
     1             IER,IWM(151))
        END IF
        NRE=NRE+1
        IF (IRES.GT.1) GOTO 600
        DO I=1,N
          DO K=1,LENP
            IF((IWM(150+K).EQ.I).AND.(IWM(150+LENP+K).EQ.J)) THEN
              WM(2+K)=(RTEM(I)-SAVR(I))*FAC
              GOTO 10
            END IF
          END DO
 10     END DO
        Y(J)=YJ
      END DO
      DO J=1,N
        SJ=S(J)
        R=DMAX1(SRUR*DABS(SJ),0.01D0/EWT(J))
        S(J)=S(J)+R
        FAC=-1.D0/R
        IF (MITER.NE.6) THEN
          CALL RES(NEQ,T,Y,S,RTEM,RPAR,LRP,IPAR,LIP)
        ELSE
          CALL RES(NEQ,T,Y,S,RTEM,WM(3),RPAR,LRP,IPAR,LIP,0,CON,
     1             IER,IWM(151))
        END IF
        NRE=NRE+1
        IF (IRES.GT.1) GOTO 600
        DO I=1,N
          DO K=1,LENP
            IF((IWM(150+K).EQ.I).AND.(IWM(150+LENP+K).EQ.J)) THEN
              WM(2+K)=WM(2+K)+(RTEM(I)-SAVR(I))*FAC
              GOTO 20
            END IF
          END DO
 20     END DO
        S(J)=SJ
      END DO
      IRES=1
      IF (MITER.NE.6) THEN
        CALL RES(NEQ,T,Y,S,SAVR,RPAR,LRP,IPAR,LIP)
      ELSE
        CALL RES(NEQ,T,Y,S,SAVR,WM(3),RPAR,LRP,IPAR,LIP,0,CON,IER,
     1           IWM(151))
      END IF
      NRE=NRE+1
 600  RETURN
      END
      SUBROUTINE OPPRINT(T,INDEX,LENP,P,CHOUT)
      IMPLICIT NONE
      DOUBLE PRECISION T,P(1)
      INTEGER CHOUT,INDEX(1),I,LENP

      WRITE(CHOUT,100) T
      WRITE(CHOUT,*) '--------------------------------',
     1               '-----'
      WRITE(CHOUT,*) '    ROW       COLUMN       VALUE'
      WRITE(CHOUT,*) '--------------------------------',
     1               '-----'
      DO I=1,LENP
        WRITE(CHOUT,110) INDEX(I),INDEX(LENP+I),
     1  P(I)
      END DO
      WRITE(CHOUT,*)
 100  FORMAT('* TIME =', G13.6)
 110  FORMAT(I8,'  ',I8,'       ',G13.6)
      RETURN
      END
      SUBROUTINE ROOTS (NG,HMIN,JFLAG,X0,X1,G0,G1,GX,X,JROOT,ALPHA,X2,
     1                  IMAX,LAST)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par : 
C 
C     Procedure appelee :
C---------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)

      INTEGER NG,JFLAG,JROOT(*),IMAX,LAST,I,IMXOLD,NXLAST
      DOUBLE PRECISION HMIN,X0,X1,G0(*),G1(*),GX(*),X,ALPHA,X2,T2,
     1                 TMAX
      LOGICAL ZROOT,SGNCHG,XROOT
                                                         
                           
C-----------------------------------------------------------------------        
C THIS SUBROUTINE FINDS THE LEFTMOST ROOT OF A SET OF ARBITRARY                 
C FUNCTIONS GI(X) (I = 1,...,NG) IN AN INTERVAL (X0,X1).  ONLY ROOTS            
C OF ODD MULTIPLICITY (I.E. CHANGES OF SIGN OF THE GI) ARE FOUND.               
C HERE THE SIGN OF X1 - X0 IS ARBITRARY, BUT IS CONSTANT FOR A GIVEN            
C PROBLEM, AND -LEFTMOST- MEANS NEAREST TO X0.                                  
C THE VALUES OF THE VECTOR-VALUED FUNCTION G(X) = (GI, I=1...NG)                
C ARE COMMUNICATED THROUGH THE CALL SEQUENCE OF ROOTS.                          
C THE METHOD USED IS THE ILLINOIS ALGORITHM.                                    
C                                                                               
C REFERENCE..                                                                   
C KATHIE L. HIEBERT AND LAWRENCE F. SHAMPINE, IMPLICITLY DEFINED                
C OUTPUT POINTS FOR SOLUTIONS OF ODE-S, SANDIA REPORT SAND80-0180,              
C FEBRUARY, 1980.                                                               
C                                                                               
C DESCRIPTION OF PARAMETERS.                                                    
C                                                                               
C NG     = NUMBER OF FUNCTIONS GI, OR THE NUMBER OF COMPONENTS OF               
C          THE VECTOR VALUED FUNCTION G(X).  INPUT ONLY.                        
C                                                                               
C HMIN   = RESOLUTION PARAMETER IN X.  INPUT ONLY.  WHEN A ROOT IS              
C          FOUND, IT IS LOCATED ONLY TO WITHIN AN ERROR OF HMIN IN X.           
C          TYPICALLY, HMIN SHOULD BE SET TO SOMETHING ON THE ORDER OF           
C               100 * UROUND * MAX(ABS(X0),ABS(X1)),                            
C          WHERE UROUND IS THE UNIT ROUNDOFF OF THE MACHINE.                    
C                                                                               
C JFLAG  = INTEGER FLAG FOR INPUT AND OUTPUT COMMUNICATION.                     
C                                                                               
C          ON INPUT, SET JFLAG = 0 ON THE FIRST CALL FOR THE PROBLEM,           
C          AND LEAVE IT UNCHANGED UNTIL THE PROBLEM IS COMPLETED.               
C          (THE PROBLEM IS COMPLETED WHEN JFLAG .GE. 2 ON RETURN.)              
C                                                                               
C          ON OUTPUT, JFLAG HAS THE FOLLOWING VALUES AND MEANINGS..             
C          JFLAG = 1 MEANS ROOTS NEEDS A VALUE OF G(X).  SET GX = G(X)          
C                    AND CALL ROOTS AGAIN.                                      
C          JFLAG = 2 MEANS A ROOT HAS BEEN FOUND.  THE ROOT IS                  
C                    AT X, AND GX CONTAINS G(X).  (ACTUALLY, X IS THE           
C                    RIGHTMOST APPROXIMATION TO THE ROOT ON AN INTERVAL         
C                    (X0,X1) OF SIZE HMIN OR LESS.)                             
C          JFLAG = 3 MEANS X = X1 IS A ROOT, WITH ONE OR MORE OF THE GI         
C                    BEING ZERO AT X1 AND NO SIGN CHANGES IN (X0,X1).           
C                    GX CONTAINS G(X) ON OUTPUT.                                
C          JFLAG = 4 MEANS NO ROOTS (OF ODD MULTIPLICITY) WERE                  
C                    FOUND IN (X0,X1) (NO SIGN CHANGES).                        
C                                                                               
C X0,X1  = ENDPOINTS OF THE INTERVAL WHERE ROOTS ARE SOUGHT.                    
C          X1 AND X0 ARE INPUT WHEN JFLAG = 0 (FIRST CALL), AND                 
C          MUST BE LEFT UNCHANGED BETWEEN CALLS UNTIL THE PROBLEM IS            
C          COMPLETED.  X0 AND X1 MUST BE DISTINCT, BUT X1 - X0 MAY BE           
C          OF EITHER SIGN.  HOWEVER, THE NOTION OF -LEFT- AND -RIGHT-           
C          WILL BE USED TO MEAN NEARER TO X0 OR X1, RESPECTIVELY.               
C          WHEN JFLAG .GE. 2 ON RETURN, X0 AND X1 ARE OUTPUT, AND               
C          ARE THE ENDPOINTS OF THE RELEVANT INTERVAL.                          
C                                                                               
C G0,G1  = ARRAYS OF LENGTH NG CONTAINING THE VECTORS G(X0) AND G(X1),          
C          RESPECTIVELY.  WHEN JFLAG = 0, G0 AND G1 ARE INPUT AND               
C          NONE OF THE G0(I) SHOULD BE BE ZERO.                                 
C          WHEN JFLAG .GE. 2 ON RETURN, G0 AND G1 ARE OUTPUT.                   
C                                                                               
C GX     = ARRAY OF LENGTH NG CONTAINING G(X).  GX IS INPUT                     
C          WHEN JFLAG = 1, AND OUTPUT WHEN JFLAG .GE. 2.                        
C                                                                               
C X      = INDEPENDENT VARIABLE VALUE.  OUTPUT ONLY.                            
C          WHEN JFLAG = 1 ON OUTPUT, X IS THE POINT AT WHICH G(X)               
C          IS TO BE EVALUATED AND LOADED INTO GX.                               
C          WHEN JFLAG = 2 OR 3, X IS THE ROOT.                                  
C          WHEN JFLAG = 4, X IS THE RIGHT ENDPOINT OF THE INTERVAL, X1.         
C                                                                               
C JROOT  = INTEGER ARRAY OF LENGTH NG.  OUTPUT ONLY.                            
C          WHEN JFLAG = 2 OR 3, JROOT INDICATES WHICH COMPONENTS                
C          OF G(X) HAVE A ROOT AT X.  JROOT(I) IS 1 IF THE I-TH                 
C          COMPONENT HAS A ROOT, AND JROOT(I) = 0 OTHERWISE.                    
C                                                                               
C NOTE.. THIS ROUTINE USES THE COMMON BLOCK /LSR001/ TO SAVE                    
C THE VALUES OF CERTAIN VARIABLES BETWEEN CALLS (OWN VARIABLES).                
C-----------------------------------------------------------------------        

C JFLAG .NE. 1.  CHECK FOR CHANGE IN SIGN OF G OR ZERO AT X1. ----------

      IF (JFLAG.EQ.1) GOTO 200
      IMAX=0
      TMAX=ZERO
      ZROOT=.FALSE.
      DO I=1,NG
        IF (DABS(G1(I)).LE.ZERO) THEN
          ZROOT=.TRUE.
          ELSE

C AT THIS POINT, G0(I) HAS BEEN CHECKED AND CANNOT BE ZERO. ------------

          IF (DSIGN(1.0D0,G0(I)).NE.DSIGN(1.0D0,G1(I))) THEN
            T2=DABS(G1(I)/(G1(I)-G0(I)))
            IF (T2.GT.TMAX) THEN
              TMAX=T2
              IMAX=I
            END IF
          END IF
        END IF
      END DO
      SGNCHG=.FALSE.
      IF (IMAX.GT.0) SGNCHG=.TRUE.
      IF (.NOT.SGNCHG) GOTO 400

C THERE IS A SIGN CHANGE.  FIND THE FIRST ROOT IN THE INTERVAL. --------

      XROOT=.FALSE.
      NXLAST=0
      LAST=1

C REPEAT UNTIL THE FIRST ROOT IN THE INTERVAL IS FOUND.  LOOP POINT. ---

 150  CONTINUE
      IF (XROOT) GOTO 300
      IF (NXLAST.NE.LAST) THEN
        ALPHA=1.0D0
        ELSE
        IF (LAST.EQ.0) THEN
          ALPHA=2.0D0*ALPHA
          ELSE
          ALPHA=0.5D0*ALPHA
        END IF
      END IF
      X2=X1-(X1-X0)*G1(IMAX)/(G1(IMAX)-ALPHA*G0(IMAX))
      IF ((DABS(X2-X0).LT.HMIN).AND.
     1   (DABS(X1-X0).GT.10.0D0*HMIN)) X2=X0+0.1D0*(X1-X0)
      JFLAG=1
      X=X2

C RETURN TO THE CALLING ROUTINE TO GET A VALUE OF GX = G(X). -----------

      RETURN

C CHECK TO SEE IN WHICH INTERVAL G CHANGES SIGN. -----------------------

 200  IMXOLD=IMAX
      IMAX=0
      TMAX=ZERO
      ZROOT=.FALSE.
      DO I=1,NG
        IF (DABS(GX(I)).LE.ZERO) THEN
          ZROOT=.TRUE.
          ELSE

C NEITHER G0(I) NOR GX(I) CAN BE ZERO AT THIS POINT. -------------------

          IF (DSIGN(1.0D0,G0(I)).NE.DSIGN(1.0D0,GX(I))) THEN
            T2=DABS(GX(I)/(GX(I)-G0(I)))
            IF (T2.GT.TMAX) THEN
              TMAX=T2
              IMAX=I
            END IF
          END IF
        END IF
      END DO
      IF (IMAX.GT.0) THEN
        SGNCHG=.TRUE.
        ELSE
        SGNCHG=.FALSE.
        IMAX=IMXOLD
      END IF
      NXLAST = LAST
      IF (SGNCHG) THEN

C SIGN CHANGE BETWEEN X0 AND X2, SO REPLACE X1 WITH X2. ----------------

        X1=X2
        CALL DCOPY(NG,GX,1,G1,1)
        LAST=1
        XROOT=.FALSE.
        ELSE
        IF (ZROOT) THEN

C ZERO VALUE AT X2 AND NO SIGN CHANGE IN (X0,X2), SO X2 IS A ROOT. -----

          X1=X2
          CALL DCOPY(NG,GX,1,G1,1)
          XROOT=.TRUE.
          ELSE

C NO SIGN CHANGE BETWEEN X0 AND X2.  REPLACE X0 WITH X2. ---------------

          CALL DCOPY(NG,GX,1,G0,1)
          X0=X2
          LAST=0
          XROOT=.FALSE.
        END IF
      END IF
      IF (DABS(X1-X0).LE.HMIN) XROOT=.TRUE.
      GOTO 150

C RETURN WITH X1 AS THE ROOT.  SET JROOT.  SET X = X1 AND GX = G1. -----

 300  JFLAG=2
      X=X1
      CALL DCOPY(NG,G1,1,GX,1)
      DO I=1,NG
        JROOT(I)=0
        IF ((DABS(G1(I)).LE.ZERO).OR.
     1     (DSIGN(1.0D0,G0(I)).NE.DSIGN(1.0D0,G1(I)))) JROOT(I)=1
      END DO
      RETURN

C NO SIGN CHANGE IN THE INTERVAL.  CHECK FOR ZERO AT RIGHT ENDPOINT. ---

 400  IF (ZROOT) THEN

C ZERO VALUE AT X1 AND NO SIGN CHANGE IN (X0,X1).  RETURN JFLAG = 3. ---

        X=X1
        CALL DCOPY(NG,G1,1,GX,1)
        DO I=1,NG
          JROOT(I)=0
          IF (DABS(G1(I)).LE.ZERO) JROOT(I)=1
        END DO
        JFLAG=3
        ELSE

C NO SIGN CHANGES IN THIS INTERVAL.  SET X = X1, RETURN JFLAG = 4. -----

        CALL DCOPY(NG,G1,1,GX,1)
        X=X1
        JFLAG=4
      END IF
      RETURN
C----------------------- Fin de la routine ROOTS -----------------------        
      END
      SUBROUTINE SAINVG (T,Y,YDOTI,NEQ,RES,ADDA,JAC,H,EWT,IDID,MITER,
     1   SAVF,YH,NYH,DELTA,WM,IWM,HMIN,UROUND,RPAR,LRP,IPAR,LIP,NRE,
     2   NJE,CNTL,RINFO,KEEP,INFO,ICNTL)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C---------------------------------------------------------------------
C     SAINVG EFFECTUE UN PAS DE TAILLE H OU PLUS PETIT
C     AVEC LA METHODE D'EULER RETROGRADE, POUR
C     TROUVER YDOTI AU TEMPS INITIAL T. UNE METHODE
C     DE NEWTON MODIFIEE AVEC RELAXATION EST UTILISEE
C     POUR RESOUDRE LE CORRECTEUR.
C
C     THE PARAMETERS REPRESENT:
c     T --         independent variable
c     y --         solution vector at T
c     YDOTI --     derivative of solution vector
c     NEQ --       number of equations
c     H --         stepsize. imder may use a stepsize
c                  smaller than h.
c     EWT --       vector of weights for error
c                  criterion
c     IDID --      completion code with the following meanings
c                  idid= 1 -- yprime was found successfully
c                  idid=-12 -- ddaini failed to find yprime
c     RPAR --      real parameter array
c     IPAR --      integer parameter array
c                  that is not altered by ddaini
c     YH --        work space for ddaini
c     DELTA --     work space for ddaini
c     WM,IWM --    real and integer arrays storing
c                  matrix information
c
c     SUBROUTINE APPELEES:
c     RES, JAC, ADDA, OPNUM, UMD2FA, UMD2RF, SSOLSY 
c     ICOPY, OPPRINT, DSCAL, DAXPY, DCOPY
c
c     DERNIERES MODIFICATIONS:   le 30/07/97 par FOED NOUR
C---------------------------------------------------------------------


      IMPLICIT NONE
    
      INTEGER IWM(*),NEQ(*),KEEP(*),INFO(*),ICNTL(*),MAXIT,MJAC,N,
     1        IDID,NCF,NSF,IPAS,JCALC,M,NRE,IRES,LENP,MITER,IER,NJE,I,
     2        NYH,IERSL,IPAR(*),LRP,LIP,IONE,IZERO
      DOUBLE PRECISION Y(*),YDOTI(*),EWT(*),YH(NYH,*),DELTA(*),SAVF(*),
     1        WM(*),RPAR(*),CNTL(*),RINFO(*),DAMP,TN,T,H,YNORM,CON,S,
     2        DELNRM,UROUND,OLDNRM,RATE,HMIN,VNORM,TCONV,TRATE,ZERO,
     3        MONE,ONE
      LOGICAL CONVGD,CORRECTOR

      EXTERNAL RES,ADDA,JAC

c***************************************************
c*                    Block 1.                     *
c*               Initializations.                  *
c***************************************************

	IZERO=0
	IONE=1
	ZERO=0.0D0
	MONE=-1.D0
	ONE=1.D0

      MAXIT=25
      MJAC=1
      TCONV=0.05D0
      TRATE=0.90D0
      DAMP=1.0D0
      N=NEQ(1)
      IDID=0
      NCF=0
      NSF=0
      IPAS=1
      YNORM = VNORM(N, Y, EWT)
      CALL DCOPY(N,Y,IONE,YH(1,1),IONE)
      CALL DCOPY(N,YDOTI,IONE,YH(1,2),IONE)

c**************************************************
c*                    Block 2.                    * 
c*           Do one backward euler step.          *
c**************************************************
      DO WHILE(IPAS.GE.0)
        TN=T+H       
        CALL DAXPY(N,H,YDOTI,IONE,Y,IONE)
        JCALC=-1 
        M=0
        CONVGD=.TRUE.
        CORRECTOR=.FALSE.
        CALL DSCAL(N,ZERO,SAVF,IONE)

CFN********* DEBUT DE BOUCLE DE CONVERGENCE DU CORRECTEUR ******************

        DO WHILE((.NOT.CORRECTOR).AND.CONVGD)
          IRES=0
          CON=-H
          S=1000000.D0
          LENP=IWM(4)
          CALL DSCAL(2*LENP,ZERO,WM(3),IONE)

CFN********* CALCUL DES RESIDUS ET DU JACOBIEN *******************************

          IF(JCALC.NE.-1) THEN
            IF(MITER.EQ.5.OR.MITER.EQ.6) THEN
              CALL RES(NEQ,TN,Y,YDOTI,DELTA,WM(3),RPAR,LRP,IPAR,LIP,
     1            IZERO,CON,IER,IWM(151))
              NRE=NRE+1       
            ELSE
              CALL RES(NEQ,TN,Y,YDOTI,DELTA,RPAR,LRP,IPAR,LIP)
              NRE=NRE+1       
            END IF
          ELSE IF(MITER.EQ.1.OR.MITER.EQ.3.OR.MITER.EQ.7) THEN
            CALL RES(NEQ,TN,Y,YDOTI,DELTA,RPAR,LRP,IPAR,LIP)
            NRE=NRE+1       
            CALL JAC(NEQ,TN,Y,YDOTI,WM(3),IWM(3),IWM(1),IWM(2),RPAR,
     1                LRP,IPAR,LIP)
            NJE=NJE+1
            CALL DSCAL(LENP,CON,WM(3),IONE)
            CALL ADDA(NEQ,TN,Y,YDOTI,WM(3+LENP),IWM(3),IWM(1),IWM(2),
     1                RPAR,LRP,IPAR,LIP)
            CALL DAXPY(LENP,MONE,WM(3+LENP),IONE,WM(3),IONE)
          ELSE IF(MITER.EQ.5) THEN
            CALL RES(NEQ,TN,Y,YDOTI,DELTA,WM(3),RPAR,LRP,IPAR,LIP,IONE,
     1              CON,IER,IWM(151))
            NRE=NRE+1       
            NJE=NJE+1
          ELSE
            CALL OPNUM(NEQ,Y,TN,CON,EWT,SAVF,DELTA,YDOTI,WM,IWM,
     1               RES,RPAR,LRP,IPAR,LIP,NRE,MITER)
          END IF 

CFN********* DECOMPOSITION LU *******************************

          IF (IRES.GE.0) THEN
            IF (IWM(9).NE.0) CALL OPPRINT(TN,IWM(151),IWM(4),WM(3),
     1                                    IWM(9))
            CALL ICOPY(2*IWM(4),IWM(151),IONE,IWM(151+2*IWM(4)),IONE)
            CALL UMD2RF(N,IWM(4),IONE,.FALSE.,IWM(20),IWM(19),
     1        WM(3),IWM(151+2*IWM(4)),KEEP,CNTL,ICNTL,INFO,RINFO)
            IF(INFO(1).LT.0) THEN
              CALL ICOPY(2*IWM(4),IWM(151),IONE,IWM(151+2*IWM(4)),IONE)
              CALL UMD2FA(N,IWM(4),IONE,.FALSE.,IWM(20),IWM(19),
     1          WM(3),IWM(151+2*IWM(4)),KEEP,CNTL,ICNTL,INFO,RINFO)
            END IF
            IER=INFO(1) 
            IF ((IER.LT.0).OR.(IER.GE.4)) CONVGD=.FALSE.
          ELSE
            CONVGD=.FALSE.
          END IF

CFN ************* RESOLUTION A*X=B PAR SSOLSY ***************

        IF(CONVGD) THEN
          CALL DSCAL(N,DAMP,DELTA,IONE)
          CALL SSOLSY(NEQ,WM,IWM,DELTA,CNTL,
     1                RINFO,KEEP,INFO,ICNTL,IERSL)
          DELNRM=VNORM(N,DELTA,EWT)*DABS(H)
          CALL DAXPY(N,ONE,DELTA,IONE,YDOTI,IONE)
          CALL DSCAL(N,H,DELTA,IONE)
          CALL DAXPY(N,ONE,DELTA,IONE,Y,IONE)
       
CFN ************ TEST SUR LA VITESSE DE CONVERGENCE ********

          IF (DELNRM.LE.100.D0*UROUND*YNORM) THEN
            CORRECTOR=.TRUE.
          ELSE 
            IF (M.EQ.0) THEN
              OLDNRM=DELNRM
            ELSE
              RATE=(DELNRM/OLDNRM)**(1.0D0/DFLOAT(M))
              IF (RATE.GT.TRATE.AND.M.LT.3) RATE=TRATE
              IF (RATE.GT.TRATE) THEN
                CONVGD=.FALSE.
              ELSE
                S=RATE/(1.0D0-RATE)
              END IF
            END IF
            IF (CONVGD) THEN
              IF (S*DELNRM .LE. TCONV) THEN
                CORRECTOR=.TRUE.
              ELSE
                M=M+1
                IF (M.GE.MAXIT) THEN
                  CONVGD=.FALSE.
                ELSE IF ((M/MJAC)*MJAC.EQ.M) THEN
                  JCALC=-1 
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO

c*****************************************************
c*                   Block 3.                        *
c*        The corrector iteration converged.         *
c*                Do error test.                     *
c*****************************************************

       IF (CONVGD) THEN
         DO I=1,N
            DELTA(I)=Y(I)-YH(I,1)+(YH(I,2)-YDOTI(I))*H
         END DO
         IF(IPAS.EQ.0) THEN
           CALL DCOPY(N,YH(1,1),IONE,Y,IONE) 
         ELSE
           T=TN
           CALL DCOPY(N,Y,IONE,YH(1,1),IONE)  
           CALL DCOPY(N,YH(1,2),IONE,YDOTI,IONE)  
         END IF
         IPAS=IPAS-1
     
c****************************************************
c*                    Block 4.                       *
c*        The backward euler step failed.            *  
c*     Restore and YDOTI to their original values.   *
c*     Reduce stepsize and try again, if possible.   *
c*****************************************************

       ELSE
         CALL DCOPY(N,YH(1,1),IONE,Y,IONE)  
         CALL DCOPY(N,YH(1,2),IONE,YDOTI,IONE)  
         IF (IER.LT.0) THEN
           NSF=NSF+1
           H=H*0.25D0
           IF (NSF.GE.6.OR.DABS(H).LT.HMIN) IDID=-10
         ELSE 
           IF (IRES.GT.-2) THEN
             NCF=NCF+1
             H=H*0.25D0
             IF (NCF.GE.10.OR.DABS(H).LT.HMIN) IDID=-12
           ELSE
             IDID=-11
           END IF
         END IF
         IF(IDID.LE.-10) RETURN
       END IF
      END DO
      RETURN
C-------------------- Fin de la routine SAINVG -------------------------
      END
C-------------------------------------------------------------------
       SUBROUTINE SPAIDX(NEQ,NE,LINDEX,LVALUE)
       
        IMPLICIT NONE
        INTEGER NEQ(*),NE,LINDEX,LVALUE

c       LINDEX=3*NE+36*NEQ(1)+1+100
c       LVALUE=4*NE
        LINDEX=2*NE
        LVALUE=NE
        RETURN
        END
      SUBROUTINE SPREPJI(NEQ,Y,TN,EWT,RTEM,SAVR,S,WM,IWM,EL0,H,
     1                   JCUR,MITER,RES,JAC,ADDA,RPAR,LRP,IPAR,LIP,NRE,
     2                   NJE,CNTL,RINFO,KEEP,INFO,ICNTL,IERPJ)
C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 24 Juillet 1997
C     Modif'      : A.S. 14 Octobre 1997
C                   Elimination des COMMONs
C     Objet	  : Calcul et Factorisation LU de l'Operateur 
C	            Dynamique
C
C     Appelee par : SSTODI
C 
C     Procedures appelees : RES
C                           JAC
C                           ADDA
C                           UMD2FA
C                           UMD2RF
C---------------------------------------------------------------------
C
C     Arguments :
C
C      NEQ()       : Nombre d'Equation                       (entree)
C      Y(NEQ)      : Variables Predites                      (entree)
C      YH(NYH,NEQ) : Vecteurs de Nordiesk des Variables      (entree)
C      EWT(NEQ)    : Vecteur des Erreurs                     (entree)
C      RTEM(NEQ)   : Vecteur de Travail                      (entree)
C      SAVR(NEQ)   : Evaluation des Residus                  (sortie)
C      S(NEQ)      : Derivees                                (entree)
C      WM()        : Vecteurs de Travail contenant la     (entree/sortie)
C                    decomposition LU de l'Operateur Dynamique
C      IWM()       : Vecteur de Travail                   (entree/sortie)
C      RPAR()      : Vecteurs des Parametres Utilisateurs    (entree)
C      IPAR()      : Vecteurs des Parametres Utilisateurs    (entree)
C      EL0         : =EL(1)                                  (entree)
C      IERPJ       : Drapeau d'Erreur                        (sortie)
C                    0 = Pas d'Erreurs
C                    1 = Matrice Singuliere
C      JCUR        : Drapeau Indiquant que l'Operateur       (sortie)
C                    Dynamique a ete reactualise (JCUR = 1)
C---------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      EXTERNAL RES, JAC, ADDA

      INTEGER NEQ(*),IWM(*),IERPJ,JCUR,MITER,NRE,NJE,IER,
     1        LENP,KEEP(*),INFO(*),ICNTL(*),IPAR(*),LRP,LIP

      DOUBLE PRECISION RPAR(*),Y(*),EWT(*),RTEM(*),SAVR(*),
     1                 S(*),WM(*),CNTL(*),RINFO(*),EL0,H,TN,HL0

      HL0=H*EL0
      IERPJ=0
      JCUR=1
      LENP=IWM(4)
      CALL DSCAL(2*LENP,0.0D0,WM(3),1)

      IF (MITER.EQ.1.OR.MITER.EQ.3.OR.MITER.EQ.7) THEN
        CALL RES(NEQ,TN,Y,S,SAVR,RPAR,LRP,IPAR,LIP)
        NRE=NRE+1
        CALL JAC(NEQ,TN,Y,S,WM(3),IWM(3),IWM(1),IWM(2),RPAR,LRP,IPAR,
     1           LIP)
        NJE=NJE+1
        CALL DSCAL(LENP,-HL0,WM(3),1)
        CALL ADDA(NEQ,TN,Y,S,WM(3+LENP),IWM(3),IWM(1),IWM(2),RPAR,LRP,
     1            IPAR,LIP)
        CALL DAXPY(LENP,-1.0D0,WM(3+LENP),1,WM(3),1)
      ELSE IF (MITER.EQ.5) THEN
        CALL RES(NEQ,TN,Y,S,SAVR,WM(3),RPAR,LRP,IPAR,LIP,1,-HL0,
     1           IER,IWM(151))
        NRE=NRE+1
        NJE=NJE+1
      ELSE
        CALL OPNUM(NEQ,Y,TN,-HL0,EWT,RTEM,SAVR,S,WM,IWM,RES,RPAR,LRP,
     1             IPAR,LIP,NRE,MITER)
      END IF

      IF (IWM(9).NE.0) CALL OPPRINT(TN,IWM(151),IWM(4),WM(3),IWM(9))

      CALL ICOPY(2*IWM(4),IWM(151),1,IWM(151+2*IWM(4)),1)
      CALL UMD2RF(NEQ(1),IWM(4),1,.FALSE.,IWM(20),IWM(19),WM(3),
     1            IWM(151+2*IWM(4)),KEEP,CNTL,ICNTL,INFO,RINFO)

      IF (INFO(1).LT.0) THEN
        CALL ICOPY(2*IWM(4),IWM(151),1,IWM(151+2*IWM(4)),1)
        CALL UMD2FA(NEQ(1),IWM(4),1,.FALSE.,IWM(20),IWM(19),WM(3),
     1              IWM(151+2*IWM(4)),KEEP,CNTL,ICNTL,INFO,RINFO)
      END IF

      IF ((INFO(1).LT.0).OR.(INFO(1).GE.4)) IERPJ=-1
      RETURN
C----------------------- Fin de la Routine SPREPJI ---------------------
      END
      SUBROUTINE SRCHEK(JOB,GEX,NEQ,N,Y,YH,NYH,G0,G1,GX,JROOT,IRT,YDOT,
     1                  TN,H,HU,UROUND,T0,TLAST,TOUTC,IRFND,ITASKC,L,NQ,
     2                  NGC,NGE,RPAR,LRP,IPAR,LIP,ALPHA,X2,IMAX,LAST,
     3                  MESFLG,LUNIT)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par : 
C 
C     Procedure appelee : 
C---------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      INTEGER IZERO
      PARAMETER (IZERO=0)

      EXTERNAL GEX

      INTEGER JOB,NEQ(*),NYH,JROOT(*),IRT,IRFND,ITASKC,NGC,NGE,I,IFLAG,
     1        JFLAG,L,NQ,N,IMAX,LAST,MESFLG,LUNIT,IPAR(*),LRP,LIP
      DOUBLE PRECISION Y(*),YH(NYH,*),G0(*),G1(*),GX(*),RPAR(*),YDOT(*),
     1        TN,UROUND,H,T0,TLAST,TOUTC,HMING,T1,TEMP1,TEMP2,X,HU,X2,
     2        ALPHA
      LOGICAL ZROOT

C-----------------------------------------------------------------------
C THIS ROUTINE CHECKS FOR THE PRESENCE OF A ROOT IN THE
C VICINITY OF THE CURRENT T, IN A MANNER DEPENDING ON THE
C INPUT FLAG JOB.  IT CALLS SUBROUTINE ROOTS TO LOCATE THE ROOT
C AS PRECISELY AS POSSIBLE.
C
C IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, RCHEK
C USES THE FOLLOWING FOR COMMUNICATION..
C JOB    = INTEGER FLAG INDICATING TYPE OF CALL..
C          JOB = 1 MEANS THE PROBLEM IS BEING INITIALIZED, AND RCHEK
C                  IS TO LOOK FOR A ROOT AT OR VERY NEAR THE INITIAL T.
C          JOB = 2 MEANS A CONTINUATION CALL TO THE SOLVER WAS JUST
C                  MADE, AND RCHEK IS TO CHECK FOR A ROOT IN THE
C                  RELEVANT PART OF THE STEP LAST TAKEN.
C          JOB = 3 MEANS A SUCCESSFUL STEP WAS JUST TAKEN, AND RCHEK
C                  IS TO LOOK FOR A ROOT IN THE INTERVAL OF THE STEP.
C G0     = ARRAY OF LENGTH NG, CONTAINING THE VALUE OF G AT T = T0.
C          G0 IS INPUT FOR JOB .GE. 2 AND ON OUTPUT IN ALL CASES.
C G1,GX  = ARRAYS OF LENGTH NG FOR WORK SPACE.
C IRT    = COMPLETION FLAG..
C          IRT = 0  MEANS NO ROOT WAS FOUND.
C          IRT = -1 MEANS JOB = 1 AND A ROOT WAS FOUND TOO NEAR TO T.
C          IRT = 1  MEANS A LEGITIMATE ROOT WAS FOUND (JOB = 2 OR 3).
C                   ON RETURN, T0 IS THE ROOT LOCATION, AND Y IS THE
C                   CORRESPONDING SOLUTION VECTOR.
C T0     = VALUE OF T AT ONE ENDPOINT OF INTERVAL OF INTEREST.  ONLY
C          ROOTS BEYOND T0 IN THE DIRECTION OF INTEGRATION ARE SOUGHT.
C          T0 IS INPUT IF JOB .GE. 2, AND OUTPUT IN ALL CASES.
C          T0 IS UPDATED BY RCHEK, WHETHER A ROOT IS FOUND OR NOT.
C TLAST  = LAST VALUE OF T RETURNED BY THE SOLVER (INPUT ONLY).
C TOUTC  = COPY OF TOUT (INPUT ONLY).
C IRFND  = INPUT FLAG SHOWING WHETHER THE LAST STEP TAKEN HAD A ROOT.
C          IRFND = 1 IF IT DID, = 0 IF NOT.
C ITASKC = COPY OF ITASK (INPUT ONLY).
C NGC    = COPY OF NG (INPUT ONLY).
C-----------------------------------------------------------------------
C
      IRT=0
      DO I=1,NGC
	  JROOT(I)=0
	END DO
		
      HMING=(DABS(TN)+DABS(H))*UROUND*100.0D0
 
      GOTO (100,200,300),JOB
 
C EVALUATE G AT INITIAL T, AND CHECK FOR ZERO VALUES. ------------------

 100  CONTINUE
      T0=TN
      CALL GEX(NEQ,T0,Y,YDOT,G0,NGC,RPAR,LRP,IPAR,LIP)
      NGE=1
      ZROOT=.FALSE.
      DO I=1,NGC
        IF (DABS(G0(I)).LE.0.0D0) ZROOT=.TRUE.
      END DO
      IF (.NOT.ZROOT) RETURN

C G HAS A ZERO AT T.  LOOK AT G AT T + (SMALL INCREMENT). --------------

      TEMP1=DSIGN(HMING,H)
      T0=T0+TEMP1
      TEMP2=TEMP1/H
      CALL DAXPY(N,TEMP2,YH(1,2),1,Y,1)
      CALL GEX(NEQ,T0,Y,YDOT,G0,NGC,RPAR,LRP,IPAR,LIP)
      NGE=NGE+1
      ZROOT=.FALSE.
      DO I=1,NGC
        IF (DABS(G0(I)).LE.0.0D0) ZROOT=.TRUE.
      END DO
      IF (.NOT.ZROOT) RETURN

C G HAS A ZERO AT T AND ALSO CLOSE TO T.  TAKE ERROR RETURN. -----------

      IRT=-1
      RETURN
 
 200  CONTINUE
      IF (IRFND.EQ.0) GOTO 260

C IF A ROOT WAS FOUND ON THE PREVIOUS STEP, EVALUATE G0 = G(T0). -------

      CALL INTDY(N,T0,TN,0,YH,NYH,Y,H,HU,UROUND,L,NQ,IFLAG,MESFLG,
     1           LUNIT)
      CALL GEX(NEQ,T0,Y,YDOT,G0,NGC,RPAR,LRP,IPAR,LIP)
      NGE=NGE+1
      ZROOT=.FALSE.
      DO I=1,NGC
        IF (DABS(G0(I)).LE.0.0D0) ZROOT=.TRUE.
      END DO
      IF (.NOT.ZROOT) GOTO 260

C G HAS A ZERO AT T0.  LOOK AT G AT T + (SMALL INCREMENT). -------------

      TEMP1=DSIGN(HMING,H)
      T0=T0+TEMP1
      IF ((T0-TN)*H.LT.0.0D0) 
     1  CALL INTDY(N,T0,TN,0,YH,NYH,Y,H,HU,UROUND,L,NQ,IFLAG,
     2             MESFLG,LUNIT)
      TEMP2=TEMP1/H
      CALL DAXPY(N,TEMP2,YH(1,2),1,Y,1)
      CALL GEX(NEQ,T0,Y,YDOT,G0,NGC,RPAR,LRP,IPAR,LIP)
      NGE=NGE+1
      ZROOT=.FALSE.
      DO I=1,NGC
        IF (DABS(G0(I)).LE.0.0D0) THEN
          JROOT(I)=1
          ZROOT=.TRUE.
        END IF
      END DO
      IF (.NOT.ZROOT) GOTO 260

C G HAS A ZERO AT T0 AND ALSO CLOSE TO T0.  RETURN ROOT. ---------------

      IRT=1
      RETURN

C     HERE, G0 DOES NOT HAVE A ROOT
C G0 HAS NO ZERO COMPONENTS.  PROCEED TO CHECK RELEVANT INTERVAL. ------

 260  IF (TN.EQ.TLAST) RETURN

C SET T1 TO TN OR TOUTC, WHICHEVER COMES FIRST, AND GET G AT T1. -------

 300  CONTINUE
      IF ((ITASKC.EQ.2.OR.ITASKC.EQ.3.OR.ITASKC.EQ.5)
     1   .OR.((TOUTC-TN)*H.GE.0.0D0)) THEN
        T1=TN
        CALL DCOPY(N,YH,1,Y,1)
      ELSE
        T1=TOUTC
        IF ((T1-T0)*H.LE.0.0D0) RETURN
        CALL INTDY(N,T1,TN,0,YH,NYH,Y,H,HU,UROUND,L,NQ,IFLAG,
     1             MESFLG,LUNIT)
      END IF
      CALL GEX(NEQ,T1,Y,YDOT,G1,NGC,RPAR,LRP,IPAR,LIP)
      NGE=NGE+1

C CALL ROOTS TO SEARCH FOR ROOT IN INTERVAL FROM T0 TO T1. -------------

      JFLAG=0
      CALL ROOTS(NGC,HMING,JFLAG,T0,T1,G0,G1,GX,X,JROOT,ALPHA,X2,IMAX,
     1           LAST)
      DO WHILE(JFLAG.LE.1)
        CALL INTDY(N,X,TN,0,YH,NYH,Y,H,HU,UROUND,L,NQ,IFLAG,
     1             MESFLG,LUNIT)
        CALL GEX(NEQ,X,Y,YDOT,GX,NGC,RPAR,LRP,IPAR,LIP)
        NGE=NGE+1
        CALL ROOTS(NGC,HMING,JFLAG,T0,T1,G0,G1,GX,X,JROOT,ALPHA,X2,
     1             IMAX,LAST)
      END DO
      T0=X
      CALL DCOPY(NGC,GX,1,G0,1)
      IF (JFLAG.EQ.4) RETURN

C FOUND A ROOT.  INTERPOLATE TO X AND RETURN. --------------------------

      CALL INTDY(N,X,TN,0,YH,NYH,Y,H,HU,UROUND,L,NQ,IFLAG,
     1           MESFLG,LUNIT)
      CALL INTDY(N,X,TN,1,YH,NYH,YDOT,H,HU,UROUND,L,NQ,IFLAG,
     1           MESFLG,LUNIT)
      CALL GEX(NEQ,X,Y,YDOT,GX,NGC,RPAR,LRP,IPAR,LIP)
      IRT=1
      RETURN
C----------------------- Fin de la routine RCHEK -----------------------
      END
      SUBROUTINE SSOLSY(NEQ,WM,IWM,X,CNTL,RINFO,KEEP,
     1                  INFO,ICNTL,IERSL)
C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 24 Juillet 1997
C     Objet	  : Resolution du Systeme AX = B
C     Modif'      : A.S. 14 Octobre 1997
C                   Elimination des COMMONs
C
C     Appelee par : SAINVG, SSTODI
C 
C     Procedure appelee : UMD2SO
C---------------------------------------------------------------------
C
C     Arguments :
C
C      WM()     : Vecteurs de Travail contenant la        (entree)
C                 decomposition LU de l'Operateur Dynamique
C      IWM()    : Vecteur de Travail                      (entree)
C      X(NEQ)   : Vaut B en entree et X en sortie      (entree/sortie)
C---------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      INTEGER IWM(*),NEQ(*),IERSL,IER,KEEP(*),INFO(*),
     1        ICNTL(*)
      DOUBLE PRECISION WM(*),X(*),CNTL(*),RINFO(*)

      IERSL = 0
      IER = 0
      CALL UMD2SO(NEQ(1),0,.FALSE.,IWM(20),IWM(19),WM(3),
     1            IWM(151+2*IWM(4)),KEEP,X,X,WM(3+IWM(20)),
     2            CNTL,ICNTL,INFO,RINFO)
      IF (INFO(1).LT.0) IERSL=-1
      RETURN
C----------------------- Fin de la Routine SSOLSY -----------------------
      END
      SUBROUTINE SSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,
     1                  IWM,RW,IW,EL,ELCO,TESCO,RES,ADDA,JAC,RPAR,LRP,
     2                  IPAR,LIP,ALPHA)

C---------------------------------------------------------------------
C     
C     Edition	  : Alain Sargousse 14 Octobre 1997
C     Objet	  : 
C
C     Appelee par :
C 
C     Procedure appelee :
C------------------------------------------------------------------------
C STODI PERFORMS ONE STEP OF THE INTEGRATION OF AN INITIAL VALUE
C PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS.
C NOTE.. STODI IS INDEPENDENT OF THE VALUE OF THE ITERATION METHOD
C INDICATOR MITER, AND HENCE IS INDEPENDENT
C OF THE TYPE OF CHORD METHOD USED, OR THE JACOBIAN STRUCTURE.
C COMMUNICATION WITH STODI IS DONE WITH THE FOLLOWING VARIABLES..
C
C NEQ    = INTEGER ARRAY CONTAINING PROBLEM SIZE IN NEQ(1), AND
C          PASSED AS THE NEQ ARGUMENT IN ALL CALLS TO RES.
C Y      = AN ARRAY OF LENGTH .GE. N USED AS THE Y ARGUMENT IN
C          ALL CALLS TO RES.
C NEQ    = INTEGER ARRAY CONTAINING PROBLEM SIZE IN NEQ(1), AND
C          PASSED AS THE NEQ ARGUMENT IN ALL CALLS TO RES AND G.
C YH     = AN NYH BY LMAX ARRAY CONTAINING THE DEPENDENT VARIABLES
C          AND THEIR APPROXIMATE SCALED DERIVATIVES, WHERE
C          LMAX = MAXORD + 1.  YH(I,J+1) CONTAINS THE APPROXIMATE
C          J-TH DERIVATIVE OF Y(I), SCALED BY H**J/FACTORIAL(J)
C          (J = 0,1,...,NQ).  ON ENTRY FOR THE FIRST STEP, THE FIRST
C          TWO COLUMNS OF YH MUST BE SET FROM THE INITIAL VALUES.
C NYH    = A CONSTANT INTEGER .GE. N, THE FIRST DIMENSION OF YH.
C YH1    = A ONE-DIMENSIONAL ARRAY OCCUPYING THE SAME SPACE AS YH.
C EWT    = AN ARRAY OF LENGTH N CONTAINING MULTIPLICATIVE WEIGHTS
C          FOR LOCAL ERROR MEASUREMENTS.  LOCAL ERRORS IN Y(I) ARE
C          COMPARED TO 1.0/EWT(I) IN VARIOUS ERROR TESTS.
C SAVF   = AN ARRAY OF WORKING STORAGE, OF LENGTH N. ALSO USED FOR
C          INPUT OF YH(*,MAXORD+2) WHEN JSTART = -1 AND MAXORD IS LESS
C          THAN THE CURRENT ORDER NQ.
C          SAME AS YDOTI IN LSODI.
C SAVR   = AN ARRAY OF WORKING STORAGE, OF LENGTH N.
C ACOR   = A WORK ARRAY OF LENGTH N USED FOR THE ACCUMULATED
C          CORRECTIONS. ON A SUCCESFUL RETURN, ACOR(I) CONTAINS
C          THE ESTIMATED ONE-STEP LOCAL ERROR IN Y(I).
C WM,IWM = REAL AND INTEGER WORK ARRAYS ASSOCIATED WITH MATRIX
C          OPERATIONS IN CHORD ITERATION.
C CCMAX  = MAXIMUM RELATIVE CHANGE IN H*EL0 BEFORE SPREPJI IS CALLED.
C H      = THE STEP SIZE TO BE ATTEMPTED ON THE NEXT STEP.
C          H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING THE
C          PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE, BUT ITS
C          SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.
C HMIN   = THE MINIMUM ABSOLUTE VALUE OF THE STEP SIZE H TO BE USED.
C HMXI   = INVERSE OF THE MAXIMUM ABSOLUTE VALUE OF H TO BE USED.
C          HMXI = 0.0 IS ALLOWED AND CORRESPONDS TO AN INFINITE HMAX.
C          HMIN AND HMXI MAY BE CHANGED AT ANY TIME, BUT WILL NOT
C          TAKE EFFECT UNTIL THE NEXT CHANGE OF H IS CONSIDERED.
C TN     = THE INDEPENDENT VARIABLE. TN IS UPDATED ON EACH STEP TAKEN.
C JSTART = AN INTEGER USED FOR INPUT ONLY, WITH THE FOLLOWING
C          VALUES AND MEANINGS..
C               0  PERFORM THE FIRST STEP.
C           .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST.
C              -1  TAKE THE NEXT STEP WITH A NEW VALUE OF H, MAXORD,
C                    N, METH, MITER, AND/OR MATRIX PARAMETERS.
C              -2  TAKE THE NEXT STEP WITH A NEW VALUE OF H,
C                    BUT WITH OTHER INPUTS UNCHANGED.
C          ON RETURN, JSTART IS SET TO 1 TO FACILITATE CONTINUATION.
C KFLAG  = A COMPLETION CODE WITH THE FOLLOWING MEANINGS..
C               0  THE STEP WAS SUCCESFUL.
C              -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED.
C              -2  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED.
C              -3  RES ORDERED IMMEDIATE RETURN.
C              -4  ERROR CONDITION FROM RES COULD NOT BE AVOIDED.
C              -5  FATAL ERROR IN SPREPJI OR SSOLSY.
C          A RETURN WITH KFLAG = -1, -2, OR -4 MEANS EITHER
C          ABS(H) = HMIN OR 10 CONSECUTIVE FAILURES OCCURRED.
C          ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF TN AND
C          THE YH ARRAY ARE AS OF THE BEGINNING OF THE LAST
C          STEP, AND H IS THE LAST STEP SIZE ATTEMPTED.
C MAXORD = THE MAXIMUM ORDER OF INTEGRATION METHOD TO BE ALLOWED.
C MAXCOR = THE MAXIMUM NUMBER OF CORRECTOR ITERATIONS ALLOWED.
C MSBP   = MAXIMUM NUMBER OF STEPS BETWEEN SPREPJI CALLS.
C MXNCF  = MAXIMUM NUMBER OF CONVERGENCE FAILURES ALLOWED.
C METH/MITER = THE METHOD FLAGS.  SEE DESCRIPTION IN DRIVER.
C N      = THE NUMBER OF FIRST-ORDER DIFFERENTIAL EQUATIONS.
C-----------------------------------------------------------------------

CLLL. OPTIMIZE

      IMPLICIT NONE

      EXTERNAL RES, ADDA, JAC

      INTEGER NEQ(*),NYH,IWM(*),I,I1,IREDO,IRES,IRET,J,JB,M,
     1        NCF,NEWQ,IER,IW(*),IPAR(*),LRP,LIP
      DOUBLE PRECISION RPAR(*),Y(*),YH(NYH,*),YH1(*),EWT(*),SAVF(*),
     1        SAVR(*),ACOR(*),WM(*),EL(*),ELCO(13,*),
     2        TESCO(3,*),DCON,DDN,DEL,DELP,DSM,DUP,EXDN,EXSM,
     3        EXUP,R,RH,RHDN,RHSM,RHUP,TOLD,VNORM,RW(*),ALPHA(*),VNORMM

      LOGICAL CONTINUE,CONVCOR

C------------------------------------------------------------------------
      INTEGER LTRET,LCONIT,LCRATE,LEL,LELCO,LHOLD,LRMAX,LTESCO,LCCMAX,
     1        LEL0,LH,LHMIN,LHMXI,LHU,LRC,LTN,LUROUND,LALPHA,LX2,LT0,
     2        LTLAST,LTOUTC,LCNTL,LRINFO,LMBAND,LNE,LCHOUT,LINDEX,
     3        LVALUE,LILLIN,LINIT,LLYH,LLEWT,LLACOR,LLSAVR,LLWM,LLIWM,
     4        LMXSTEP,LMXHNIL,LNHNIL,LNTREP,LNSLAST,LNYH,LIALTH,LIPUP,
     5        LLMAX,LMEO,LNQNYH,LNSLP,LICF,LIERPJ,LIERSL,LJCUR,LJSTART,
     6        LKFLAG,LL,LMETH,LMITER,LMAXORD,LMAXCOR,LMSBP,LMXNCF,LN,
     7        LNQ,LNST,LNRE,LNJE,LNQU,LLG0,LLG1,LLGX,LIMAX,LLAST,
     8        LIRFND,LITASKC,LNGC,LNGE,LMESFLG,LLUNIT,LKEEP,LINFO,
     9        LICNTL,LIDX

      
      LTRET=21
      LCONIT=LTRET+1
      LCRATE=LCONIT+1
      LEL=LCRATE+1
      LELCO=LEL+13
      LHOLD=LELCO+156
      LRMAX=LHOLD+1
      LTESCO=LRMAX+1
      LCCMAX=LTESCO+36
      LEL0=LCCMAX+1
      LH=LEL0+1
      LHMIN=LH+1
      LHMXI=LHMIN+1
      LHU=LHMXI+1
      LRC=LHU+1
      LTN=LRC+1
      LUROUND=LTN+1
      LALPHA=LUROUND+1
      LX2=LALPHA+1
      LT0=LX2+1
      LTLAST=LT0+1
      LTOUTC=LTLAST+1
      LCNTL=LTOUTC+1
      LRINFO=LCNTL+10

      LMBAND=3
      LNE=4
      LCHOUT=8
      LINDEX=19
      LVALUE=20
      LILLIN=21
      LINIT=LILLIN+1
      LLYH=LINIT+1
      LLEWT=LLYH+1
      LLACOR=LLEWT+1
      LLSAVR=LLACOR+1
      LLWM=LLSAVR+1
      LLIWM=LLWM+1
      LMXSTEP=LLIWM+1
      LMXHNIL=LMXSTEP+1
      LNHNIL=LMXHNIL+1
      LNTREP=LNHNIL+1
      LNSLAST=LNTREP+1
      LNYH=LNSLAST+1
      LIALTH=LNYH+1
      LIPUP=LIALTH+1
      LLMAX=LIPUP+1
      LMEO=LLMAX+1
      LNQNYH=LMEO+1
      LNSLP=LNQNYH+1
      LICF=LNSLP+1
      LIERPJ=LICF+1
      LIERSL=LIERPJ+1
      LJCUR=LIERSL+1
      LJSTART=LJCUR+1
      LKFLAG=LJSTART+1
      LL=LKFLAG+1
      LMETH=LL+1
      LMITER=LMETH+1
      LMAXORD=LMITER+1
      LMAXCOR=LMAXORD+1
      LMSBP=LMAXCOR+1
      LMXNCF=LMSBP+1
      LN=LMXNCF+1
      LNQ=LN+1
      LNST=LNQ+1
      LNRE=LNST+1
      LNJE=LNRE+1
      LNQU=LNJE+1
      LLG0=LNQU+1
      LLG1=LLG0+1
      LLGX=LLG1+1
      LIMAX=LLGX+1
      LLAST=LIMAX+1
      LIRFND=LLAST+1
      LITASKC=LIRFND+1
      LNGC=LITASKC+1
      LNGE=LNGC+1
      LMESFLG=LNGE+1
      LLUNIT=LMESFLG+1
      LKEEP=LLUNIT+1
      LINFO=LKEEP+20
      LICNTL=LINFO+40
      LIDX=LICNTL+20
C------------------------------------------------------------------------

      IW(LKFLAG)=0
      TOLD=RW(LTN)
      NCF=0
      IW(LIERPJ)=0
      IW(LIERSL)=0
      IW(LJCUR)=0
      IW(LICF)=0
      DELP=0.0D0
      IF (IW(LJSTART).GT.0) GOTO 200
      IF (IW(LJSTART).EQ.-1) GOTO 100
      IF (IW(LJSTART).EQ.-2) GOTO 160
C-----------------------------------------------------------------------
C ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE
C INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED
C IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL
C INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE
C OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2
C FOR THE NEXT INCREASE.
C-----------------------------------------------------------------------
      IW(LLMAX)=IW(LMAXORD)+1
      IW(LNQ)=1
      IW(LL)=2
      IW(LIALTH)=2
      RW(LRMAX)=10000.0D0
      RW(LRC)=0.0D0
      RW(LEL0)=1.0D0
      RW(LCRATE)=0.7D0
      RW(LHOLD)=RW(LH)
      IW(LMEO)=IW(LMETH)
      IW(LNSLP)=0
      IW(LIPUP)=IW(LMITER)
      IRET=3
      GOTO 140
C-----------------------------------------------------------------------
C THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1.
C IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE.
C IF AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1),
C IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP.
C IF THE CALLER HAS CHANGED METH, CFODE IS CALLED TO RESET
C THE COEFFICIENTS OF THE METHOD.
C IF THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
C ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY.
C IF H IS TO BE CHANGED, YH MUST BE RESCALED.
C IF H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1
C TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS.
C-----------------------------------------------------------------------
 100  IW(LIPUP)=IW(LMITER)
      IW(LLMAX)=IW(LMAXORD)+1
      IF (IW(LIALTH).EQ.1) IW(LIALTH)=2
      IF (IW(LMETH).EQ.IW(LMEO)) GOTO 110
      CALL CFODE(IW(LMETH),ELCO,TESCO)
      IW(LMEO)=IW(LMETH)
      IF (IW(LNQ).GT.IW(LMAXORD)) GOTO 120
      IW(LIALTH)=IW(LL)
      IRET=1
      GOTO 150
 110  IF (IW(LNQ).LE.IW(LMAXORD)) GOTO 160
 120  IW(LNQ)=IW(LMAXORD)
      IW(LL)=IW(LLMAX)
      CALL DCOPY(IW(LL),ELCO(1,IW(LNQ)),1,EL,1)
      IW(LNQNYH)=IW(LNQ)*NYH
      RW(LRC)=RW(LRC)*EL(1)/RW(LEL0)
      RW(LEL0)=EL(1)
      RW(LCONIT)=0.5D0/DFLOAT(IW(LNQ)+2)
      DDN=VNORM(IW(LN),SAVF,EWT)/TESCO(1,IW(LL))
      EXDN=1.0D0/DFLOAT(IW(LL))
      RHDN=1.0D0/(1.3D0*DDN**EXDN+0.0000013D0)
      RH=DMIN1(RHDN,1.0D0)
      IREDO=3
      IF (RW(LH).EQ.RW(LHOLD)) GOTO 170
      RH=DMIN1(RH,DABS(RW(LH)/RW(LHOLD)))
      RW(LH)=RW(LHOLD)
      GOTO 175
C-----------------------------------------------------------------------
C CFODE IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE
C CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET
C WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM.
C-----------------------------------------------------------------------
 140  CALL CFODE(IW(LMETH),ELCO,TESCO)
 150  CALL DCOPY(IW(LL),ELCO(1,IW(LNQ)),1,EL,1)
      IW(LNQNYH)=IW(LNQ)*NYH
      RW(LRC)=RW(LRC)*EL(1)/RW(LEL0)
      RW(LEL0)=EL(1)
      RW(LCONIT)=0.5D0/DFLOAT(IW(LNQ)+2)
      GOTO (160,170,200), IRET
C-----------------------------------------------------------------------
C IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
C RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO
C L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS
C FORCED BY A CONVERGENCE OR ERROR TEST FAILURE.
C-----------------------------------------------------------------------
 160  IF (RW(LH).EQ.RW(LHOLD)) GOTO 200
      RH=RW(LH)/RW(LHOLD)
      RW(LH)=RW(LHOLD)
      IREDO=3
      GOTO 175
 170  RH=DMAX1(RH,RW(LHMIN)/DABS(RW(LH)))
 175  RH=DMIN1(RH,RW(LRMAX))
      RH=RH/DMAX1(1.0D0,DABS(RW(LH))*RW(LHMXI)*RH)
      R=1.0D0
      DO J=2,IW(LL)
        R=R*RH
        CALL DSCAL(IW(LN),R,YH(1,J),1)
      END DO
      RW(LH)=RW(LH)*RH
      RW(LRC)=RW(LRC)*RH
      IW(LIALTH)=IW(LL)
      IF (IREDO.EQ.0) GOTO 690
C-----------------------------------------------------------------------
C THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY
C MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX.
C RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).
C WHEN RC DIFFERS FROM 1 BY MORE THAN CCMAX, IPUP IS SET TO MITER
C TO FORCE SPREPJI TO BE CALLED.
C IN ANY CASE, SPREPJI IS CALLED AT LEAST EVERY MSBP STEPS.
C-----------------------------------------------------------------------
 200  IF ((DABS(RW(LRC)-1.0D0).GT.RW(LCCMAX)).OR.
     1   (IW(LNST).GE.IW(LNSLP)+IW(LMSBP))) IW(LIPUP)=IW(LMITER)
      RW(LTN)=RW(LTN)+RW(LH)
      I1=IW(LNQNYH)+1
      DO JB=1,IW(LNQ)
        I1=I1-NYH
        DO I=I1,IW(LNQNYH)
          YH1(I)=YH1(I)+YH1(I+NYH)
        END DO
      END DO
C-----------------------------------------------------------------------
C UP TO MAXCOR CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS
C MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY H AND THE
C ERROR WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED
C IN ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.
C-----------------------------------------------------------------------
      CONVCOR=.TRUE.
      DO WHILE(CONVCOR)
        M=0
        CALL DAXEY(IW(LN),1.0D0/RW(LH),YH(1,2),1,SAVF,1)
        CALL DCOPY(IW(LN),YH(1,1),1,Y,1)
C-----------------------------------------------------------------------
C IF INDICATED, THE MATRIX P = A - H*EL(1)*DR/DY IS REEVALUATED AND
C PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET
C TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE.
C-----------------------------------------------------------------------
        IF (IW(LIPUP).GT.0) THEN
          IW(LIPUP)=0
          RW(LRC)=1.0D0
          RW(LCRATE)=0.7D0
          IW(LNSLP)=IW(LNST)
          CALL SPREPJI(NEQ,Y,RW(LTN),EWT,ACOR,SAVR,SAVF,
     1                 WM,IWM,RW(LEL0),RW(LH),IW(LJCUR),IW(LMITER),RES,
     2                 JAC,ADDA,RPAR,LRP,IPAR,LIP,IW(LNRE),IW(LNJE),
     3                 RW(LCNTL),RW(LRINFO),IW(LKEEP),IW(LINFO),
     4                 IW(LICNTL),IW(LIERPJ))
          IRES=IW(LIERPJ)
          IF (IW(LIERPJ).LT.0) GOTO 430
        ELSE
          IRES=1
          IF (IW(LMITER).GE.5.AND.IW(LMITER).LE.6) THEN
            CALL RES(NEQ,RW(LTN),Y,SAVF,SAVR,WM(3),RPAR,LRP,IPAR,LIP,0,
     1               0.,IER,IW(LIDX))
          ELSE
            CALL RES(NEQ,RW(LTN),Y,SAVF,SAVR,RPAR,LRP,IPAR,LIP)
          END IF
          IW(LNRE)=IW(LNRE)+1
          IRES=IABS(IRES)
        END IF
        IF (IRES.EQ.2) GOTO 435
        IF (IRES.EQ.3) GOTO 430
        CALL DSCAL(IW(LN),0.0D0,ACOR,1)
C-----------------------------------------------------------------------
C SOLVE THE LINEAR SYSTEM WITH THE CURRENT RESIDUAL AS
C RIGHT-HAND SIDE AND P AS COEFFICIENT MATRIX.
C-----------------------------------------------------------------------
        CONTINUE=.TRUE.
        DO WHILE(CONTINUE)
          CALL SSOLSY(NEQ,WM,IWM,SAVR,
     1                RW(LCNTL),RW(LRINFO),IW(LKEEP),
     2                IW(LINFO),IW(LICNTL),IW(LIERSL))
          IF (IW(LIERSL).LT.0) GOTO 430
          IF (IW(LIERSL).GT.0) THEN
            CONTINUE=.FALSE.
          ELSE
            DEL=VNORM(IW(LN),SAVR,EWT)*DABS(RW(LH))
            CALL DAXPY(IW(LN),1.0D0,SAVR,1,ACOR,1)
            DO I=1,IW(LN)
              SAVF(I)=ACOR(I)+YH(I,2)/RW(LH)
              Y(I)=YH(I,1)+EL(1)*RW(LH)*ACOR(I)
            END DO
C-----------------------------------------------------------------------
C TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE
C RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
C-----------------------------------------------------------------------
            IF (M.NE.0) RW(LCRATE)=DMAX1(0.2D0*RW(LCRATE),DEL/DELP)
            DCON=DEL*DMIN1(1.0D0,1.5D0*RW(LCRATE))/
     1           (TESCO(2,IW(LNQ))*RW(LCONIT))
            IF (DCON.LE.1.0D0) GOTO 460
            M=M+1
            IF ((M.EQ.IW(LMAXCOR)).OR.(M.GE.2.AND.DEL.GT.2.0D0*DELP)) 
     1        THEN
              CONTINUE=.FALSE.
            ELSE
              DELP=DEL
              IRES=1
              IF (IW(LMITER).GE.5) THEN
                CALL RES(NEQ,RW(LTN),Y,SAVF,SAVR,WM(3),RPAR,LRP,IPAR,
     1                   LIP,0,0.,IER,IW(LIDX))
              ELSE
                CALL RES(NEQ,RW(LTN),Y,SAVF,SAVR,RPAR,LRP,IPAR,LIP)
              ENDIF
              IW(LNRE)=IW(LNRE)+1
              IRES=IABS(IRES)
              IF (IRES.EQ.2) GOTO 435
              IF (IRES.EQ.3) CONTINUE=.FALSE.
            END IF
          END IF
        END DO
C-----------------------------------------------------------------------
C THE CORRECTORS FAILED TO CONVERGE, OR RES HAS RETURNED ABNORMALLY.
C ON A CONVERGENCE FAILURE, IF THE JACOBIAN IS OUT OF DATE, SPREPJI IS
C CALLED FOR THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS
C VALUES BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.
C TAKE AN ERROR EXIT IF IRES = 2, OR H CANNOT BE REDUCED, OR MXNCF
C FAILURES HAVE OCCURRED, OR A FATAL ERROR OCCURRED IN SPREPJI OR SSOLSY.
C-----------------------------------------------------------------------
        IW(LICF)=1
        IF (IW(LJCUR).EQ.1) THEN
          CONVCOR=.FALSE.
        ELSE
          IW(LIPUP)=IW(LMITER)
        END IF
      END DO
 430  IW(LICF)= 2
      NCF=NCF+1
      RW(LRMAX)=2.0D0
 435  RW(LTN)=TOLD
      I1=IW(LNQNYH)+1
      DO JB=1,IW(LNQ)
        I1=I1-NYH
        DO I=I1,IW(LNQNYH)
          YH1(I)=YH1(I)-YH1(I+NYH)
        END DO
      END DO
      IF (IRES.EQ.2) GOTO 680
      IF (IW(LIERPJ).LT.0.OR.IW(LIERSL).LT.0) GOTO 685
      IF (DABS(RW(LH)).LE.RW(LHMIN)*1.00001D0) GOTO 450
      IF (NCF.EQ.IW(LMXNCF)) GOTO 450
      RH=0.25D0
      IW(LIPUP)=IW(LMITER)
      IREDO=1
      GOTO 170
 450  IF (IRES.EQ.3) GOTO 680
      GOTO 670
C-----------------------------------------------------------------------
C THE CORRECTOR HAS CONVERGED.  JCUR IS SET TO 0
C TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER.
C THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500
C IF IT FAILS.
C-----------------------------------------------------------------------
 460  IW(LJCUR)=0
      IF (M.EQ.0) THEN 
        DSM=DEL/TESCO(2,IW(LNQ))
      ELSE
        DSM=DABS(RW(LH))*VNORMM(IW(LN),ACOR,EWT,ALPHA)/TESCO(2,IW(LNQ))
      END IF
      IF (DSM.GT.1.0D0) GOTO 500
C-----------------------------------------------------------------------
C AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY.
C CONSIDER CHANGING H IF IALTH = 1.  OTHERWISE DECREASE IALTH BY 1.
C IF IALTH IS THEN 1 AND NQ .LT. MAXORD, THEN ACOR IS SAVED FOR
C USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.
C IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER
C BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A
C FACTOR OF AT LEAST 1.1.  IF NOT, IALTH IS SET TO 3 TO PREVENT
C TESTING FOR THAT MANY STEPS.
C-----------------------------------------------------------------------
      IW(LKFLAG)=0
      IREDO=0
      IW(LNST)=IW(LNST)+1
      RW(LHU)=RW(LH)
      IW(LNQU)=IW(LNQ)
      DO J=1,IW(LL)
        CALL DAXPY(IW(LN),EL(J)*RW(LH),ACOR,1,YH(1,J),1)
      END DO
      IW(LIALTH)=IW(LIALTH)-1
      IF (IW(LIALTH).EQ.0) GOTO 520
      IF ((IW(LIALTH).LE.1).AND.(IW(LL).NE.IW(LLMAX)))
     1   CALL DCOPY(IW(LN),ACOR,1,YH(1,IW(LLMAX)),1)
      GOTO 700
C-----------------------------------------------------------------------
C THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.
C RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE
C TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR
C ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE
C BY A FACTOR OF 0.1 OR LESS.
C-----------------------------------------------------------------------
 500  IW(LKFLAG)=IW(LKFLAG)-1
      RW(LTN)=TOLD
      I1=IW(LNQNYH)+1
      DO JB=1,IW(LNQ)
        I1=I1-NYH
        DO I=I1,IW(LNQNYH)
          YH1(I)=YH1(I)-YH1(I+NYH)
        END DO
      END DO
      RW(LRMAX)=2.0D0
      IF (DABS(RW(LH)).LE.RW(LHMIN)*1.00001D0) GOTO 660
      IF (IW(LKFLAG).LE.-7) GOTO 660
      IREDO=2
      RHUP=0.0D0
      GOTO 540
C-----------------------------------------------------------------------
C REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS
C RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED
C AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.
C IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE.
C THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN
C ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE
C ADDITIONAL SCALED DERIVATIVE.
C-----------------------------------------------------------------------
 520  RHUP=0.0D0
      IF (IW(LL).EQ.IW(LLMAX)) GOTO 540
      DO I=1,IW(LN)
        SAVF(I)=ACOR(I)-YH(I,IW(LLMAX))
      END DO
      DUP=DABS(RW(LH))*VNORMM(IW(LN),SAVF,EWT,ALPHA)/TESCO(3,IW(LNQ))
      EXUP=1.0D0/DFLOAT(IW(LL)+1)
      RHUP=1.0D0/(1.4D0*DUP**EXUP+0.0000014D0)
 540  EXSM=1.0D0/DFLOAT(IW(LL))
      RHSM=1.0D0/(1.2D0*DSM**EXSM+0.0000012D0)
      RHDN=0.0D0
      IF (IW(LNQ).EQ.1) GOTO 560
      DDN=VNORMM(IW(LN),YH(1,IW(LL)),EWT,ALPHA)/TESCO(1,IW(LNQ))
      EXDN=1.0D0/DFLOAT(IW(LNQ))
      RHDN=1.0D0/(1.3D0*DDN**EXDN+0.0000013D0)
 560  IF (RHSM.GE.RHUP) GOTO 570
      IF (RHUP.GT.RHDN) GOTO 590
      GOTO 580
 570  IF (RHSM.LT.RHDN) GOTO 580
      NEWQ=IW(LNQ)
      RH=RHSM
      GOTO 620
 580  NEWQ=IW(LNQ)-1
      RH=RHDN
      IF (IW(LKFLAG).LT.0.AND.RH.GT.1.0D0) RH=1.0D0
      GOTO 620
 590  NEWQ=IW(LL)
      RH=RHUP
      IF (RH.LT.1.1D0) GOTO 610
      R=RW(LH)*EL(IW(LL))/DFLOAT(IW(LL))
      CALL DAXEY(IW(LN),R,ACOR,1,YH(1,NEWQ+1),1)
      GOTO 630
 610  IW(LIALTH)=3
      GOTO 700
 620  IF ((IW(LKFLAG).EQ.0).AND.(RH.LT.1.1D0)) GOTO 610
      IF (IW(LKFLAG).LE.-2) RH=DMIN1(RH,0.1D0)
C-----------------------------------------------------------------------
C IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.
C IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED.
C THEN EXIT FROM 690 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE.
C-----------------------------------------------------------------------
      IF (NEWQ.EQ.IW(LNQ)) GOTO 170
 630  IW(LNQ)=NEWQ
      IW(LL)=IW(LNQ)+1
      IRET=2
      GOTO 150
C-----------------------------------------------------------------------
C ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD
C TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.
C-----------------------------------------------------------------------
 660  IW(LKFLAG)=-1
      GOTO 720
 670  IW(LKFLAG)=-2
      GOTO 720
 680  IW(LKFLAG)=-1-IRES
      GOTO 720
 685  IW(LKFLAG)=-5
      GOTO 720
 690  RW(LRMAX)=10.0D0
 700  R=RW(LH)/TESCO(2,IW(LNQU))
      CALL DSCAL(IW(LN),R,ACOR,1)
 720  RW(LHOLD)=RW(LH)
      IW(LJSTART)=1
      RETURN
C----------------------- Fin de la Routine SSTODI ----------------------
      END
C-------------------------------------------------------------------
       SUBROUTINE TBBIDX(INDEX,NEQ,NE,LINDEX,LVALUE)
       
        IMPLICIT NONE
        INTEGER NEQ(*),INDEX(*),I,J,K,L,NE,LINDEX,LVALUE
        INTEGER LENP,ISHIFT

        LENP=NEQ(15)+NEQ(16)+NEQ(17)+NEQ(18)+NEQ(19)
        NE=LENP
c       LINDEX=3*NE+36*NEQ(1)+1+100
c       LVALUE=4*NE
        LINDEX=2*NE
        LVALUE=NE
        CALL NULLIDX(INDEX,LINDEX)
c       Creation index de E
        L=0
        DO J=NEQ(10),NEQ(12)
          DO K=1,NEQ(8)
            L=L+1
            INDEX(L)=NEQ(3)*NEQ(7)+K
            INDEX(L+LENP)=J
          END DO
        END DO
c       Creation index de A
        DO I=1,NEQ(3)
          DO J=NEQ(6),NEQ(7)
            DO K=1,NEQ(5)
              L=L+1
              IF (I.NE.1) THEN
                INDEX(L)=(I-1)*NEQ(7)+K
                INDEX(L+LENP)=(I-2)*NEQ(7)+J
              ELSE
                INDEX(L)=0
                INDEX(L+LENP)=0
              END IF
            END DO
          END DO
        END DO
c       Creation index de B
        DO I=1,NEQ(3)
          DO J=1,NEQ(7)
            DO K=1,NEQ(7)
              L=L+1
              INDEX(L)=(I-1)*NEQ(7)+K
              INDEX(L+LENP)=(I-1)*NEQ(7)+J
            END DO
          END DO
        END DO
c       Creation index de C
        DO I=1,NEQ(3)
          DO J=1,NEQ(4)
            DO K=1,NEQ(7)
              L=L+1
              IF (I.NE.NEQ(3)) THEN
                INDEX(L)=(I-1)*NEQ(7)+K
                INDEX(L+LENP)=I*NEQ(7)+J
              ELSE
                INDEX(L)=0
                INDEX(L+LENP)=0
              END IF
            END DO
          END DO
        END DO
c       Creation index de G
        DO I=1,NEQ(3)
          DO J=1,NEQ(8)
            DO K=1,NEQ(7)
              L=L+1
              INDEX(L)=(I-1)*NEQ(7)+K
              INDEX(L+LENP)=NEQ(3)*NEQ(7)+J
            END DO
          END DO
        END DO
        K=0
        DO I=1,NEQ(3)
          IF (NEQ(22+I).LT.NEQ(7)) THEN
            DO J=NEQ(22+I)+1,NEQ(7)
              DO L=1,LENP
                IF (INDEX(L).EQ.J+K)THEN
                  INDEX(L)=0
                  INDEX(L+LENP)=0
                END IF
                IF (INDEX(L+LENP).EQ.J+K) THEN
                  INDEX(L+LENP)=0
                  INDEX(L)=0
                END IF
              END DO
            END DO
          END IF
          K=K+NEQ(7)
        END DO
        K=0
        DO I=1,NEQ(3)
          ISHIFT=NEQ(7)-NEQ(22+I)
          K=K+NEQ(22+I)
          DO L=1,LENP
            IF (INDEX(L).GT.K+ISHIFT) INDEX(L)=INDEX(L)-ISHIFT
            IF (INDEX(L+LENP).GT.K+ISHIFT) INDEX(L+LENP)=
     1                             INDEX(L+LENP)-ISHIFT
          END DO
        END DO
        RETURN
        END

C---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VNORM(N,V,W)

CLLL. OPTIMIZE
C-----------------------------------------------------------------------
C THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM
C OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
C CONTAINED IN THE ARRAY W OF LENGTH N..
C   VNORM = SQRT( (1/N) * SUM( V(I)*W(I) )**2 )
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION V(*),W(*),SUM,VI,WI

      SUM=0.0D0
      DO I=1,N
        VI=DABS(V(I))
        WI=DABS(W(I))
c       IF(VI.LT.1.D-10) THEN
c          VI=1.D-10
c       ELSE IF(VI.GT.1.D+10) THEN
c          VI=1.D+10
c       END IF
c       IF(WI.LT.1.D-10) THEN
c          WI=1.D-10
c       ELSE IF(WI.GT.1.D+10) THEN
c          WI=1.D+10
c       END IF
        SUM=SUM+(VI*WI)**2
      END DO
      VNORM=DSQRT(SUM/DFLOAT(N))
      RETURN
C----------------------- Fin de la fonction VNORM ----------------------
      END

      DOUBLE PRECISION FUNCTION VNORMM(N,V,W,ALPHA)

CLLL. OPTIMIZE
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I,N1
      DOUBLE PRECISION V(*),W(*),SUM,ALPHA(*),VI,WI

      SUM=0.0D0
      N1=0
      DO I=1,N
        IF(ALPHA(I).NE.0.D0) THEN
          VI=DABS(V(I))
          WI=DABS(W(I))
c         IF(VI.LT.1.D-10) THEN
c            VI=1.D-10
c         ELSE IF(VI.GT.1.D+10) THEN
c            VI=1.D+10
c         END IF
c         IF(WI.LT.1.D-10) THEN
c            WI=1.D-10
c         ELSE IF(WI.GT.1.D+10) THEN
c            WI=1.D+10
c         END IF
          N1=N1+1
          SUM=SUM+(VI*WI*ALPHA(I))**2
        END IF
      END DO
      VNORMM=DSQRT(SUM/DFLOAT(N1))
      RETURN
C----------------------- Fin de la fonction VNORMM ---------------------
      END
      SUBROUTINE XERRWV(MSG,IERT,NI,I1,I2,NR,R1,R2,MESFLG,
     1                   LUNIT,IFIRST)

C---------------------------------------------------------------------
C
C     Edition       : Alain Sargousse 14 Octobre 1997
C     Objet   :
C
C     Appelee par :
C
C     Procedure appelee :
C---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IERT,NI,I1,I2,NR,LUN,LUNIT,MESFLG,IFIRST
      CHARACTER MSG*(*)
	DOUBLE PRECISION R1,R2
C-----------------------------------------------------------------------
C SUBROUTINES XERRWV, XSETF, AND XSETUN, AS GIVEN HERE, CONSTITUTE
C A SIMPLIFIED VERSION OF THE SLATEC ERROR HANDLING PACKAGE.
C WRITTEN BY A. C. HINDMARSH AT LLNL.  VERSION OF AUGUST 13, 1981.
C THIS VERSION IS IN DOUBLE PRECISION.
C
C ALL ARGUMENTS ARE INPUT ARGUMENTS.
C
C MSG    = THE MESSAGE (HOLLERITH LITTERAL OR INTEGER ARRAY).
C NMES   = THE LENGTH OF MSG (NUMBER OF CHARACTERS).
C NERR   = THE ERROR NUMBER (NOT USED).
C IERT   = THE ERROR TYPE..
C          1 MEANS RECOVERABLE (CONTROL RETURNS TO CALLER).
C          2 MEANS FATAL (RUN IS ABORTED--SEE NOTE BELOW).
C NI     = NUMBER OF INTEGERS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C I1,I2  = INTEGERS TO BE PRINTED, DEPENDING ON NI.
C NR     = NUMBER OF REALS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C R1,R2  = REALS TO BE PRINTED, DEPENDING ON NR.
C
C NOTE..  THIS ROUTINE IS MACHINE-DEPENDENT AND SPECIALIZED FOR USE
C IN LIMITED CONTEXT, IN THE FOLLOWING WAYS..
C 1. THE NUMBER OF HOLLERITH CHARACTERS STORED PER WORD, DENOTED
C    BY NCPW BELOW, IS A DATA-LOADED CONSTANT.
C 2. THE VALUE OF NMES IS ASSUMED TO BE AT MOST 60.
C    (MULTI-LINE MESSAGES ARE GENERATED BY REPEATED CALLS.)
C 3. IF IERT = 2, CONTROL PASSES TO THE STATEMENT   STOP
C    TO ABORT THE RUN.  THIS STATEMENT MAY BE MACHINE-DEPENDENT.
C 4. R1 AND R2 ARE ASSUMED TO BE IN DOUBLE PRECISION AND ARE PRINTED
C    IN D21.13 FORMAT.
C 5. THE COMMON BLOCK /EH0001/ BELOW IS DATA-LOADED (A MACHINE-
C    DEPENDENT FEATURE) WITH DEFAULT VALUES.
C    THIS BLOCK IS NEEDED FOR PROPER RETENTION OF PARAMETERS USED BY
C    THIS ROUTINE WHICH THE USER CAN RESET BY CALLING XSETF OR XSETUN.
C    THE VARIABLES IN THIS BLOCK ARE AS FOLLOWS..
C       MESFLG = PRINT CONTROL FLAG..
C                1 MEANS PRINT ALL MESSAGES (THE DEFAULT).
C                0 MEANS NO PRINTING.
C       LUNIT  = LOGICAL UNIT NUMBER FOR MESSAGES.
C                THE DEFAULT IS 6 (MACHINE-DEPENDENT).
C-----------------------------------------------------------------------
C THE FOLLOWING ARE INSTRUCTIONS FOR INSTALLING THIS ROUTINE
C IN DIFFERENT MACHINE ENVIRONMENTS.
C
C TO CHANGE THE DEFAULT OUTPUT UNIT, CHANGE THE DATA STATEMENT
C IN THE BLOCK DATA SUBPROGRAM BELOW.
C
C FOR A DIFFERENT NUMBER OF CHARACTERS PER WORD, CHANGE THE
C DATA STATEMENT SETTING NCPW BELOW, AND FORMAT 10.  ALTERNATIVES FOR
C VARIOUS COMPUTERS ARE SHOWN IN COMMENT CARDS.
C
C FOR A DIFFERENT RUN-ABORT COMMAND, CHANGE THE STATEMENT FOLLOWING
C STATEMENT 100 AT THE END.
C-----------------------------------------------------------------------
      IF (MESFLG.EQ.0) GOTO 100
C GET LOGICAL UNIT NUMBER. ---------------------------------------------
      LUN=LUNIT
C WRITE THE MESSAGE. ---------------------------------------------------
      IF (IFIRST.EQ.1) THEN
        WRITE(LUN,*) "---------------------------------"
        WRITE(LUN,*) " DISCo ! v1.02 beta (02/03/1998)"
        WRITE(LUN,*) "---------------------------------"
        WRITE(LUN,*)
      END IF
      WRITE (LUN,'(A)') MSG
      IF (NI.EQ.1) WRITE(LUN,20) I1
 20   FORMAT(10X,'In above message,  I1 =',I10)
      IF (NI.EQ.2) WRITE (LUN,30) I1,I2
 30   FORMAT(10X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR.EQ.1) WRITE (LUN,40) R1
 40   FORMAT(10X,'In above message,  R1 =',D21.13)
      IF (NR.EQ.2) WRITE (LUN,50) R1,R2
 50   FORMAT(10X,'In above message,  R1 =',D21.13,3X,'R2 =',D21.13)
C ABORT THE RUN IF IERT = 2. -------------------------------------------
 100  IF (IERT.NE.2) RETURN
      STOP
C----------------------- Fin de la routine XERRWV ----------------------
      END
