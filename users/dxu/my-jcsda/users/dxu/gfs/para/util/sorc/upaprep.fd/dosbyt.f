
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    DOSBYT      PACK AND WRITE A BLOCK    
C   PRGMMR: LARRY SAGER      ORG: W/NMC41    DATE: 96-08-27
C
C ABSTRACT: DOSBYT PACKS A REPORT INTO THE OUTPUT BLOCK      
C   USING SBYTES.  WHEN THE BLOCK IS FULL, IT OUTPUTS IT.  
C
C PROGRAM HISTORY LOG:
C   96-08-27  LARRY SAGER
C
C USAGE:    CALL DOSBYT(IARR,KRET,IUNO,COUT,KNEXT,IREPTS,IENT)
C   INPUT ARGUMENT LIST:
C     IARR     - REPORT                                             
C     KRET     - 12HR OLD MANDATORY LEVEL DATA
C     IUNO     - OUTPUT UNIT NUMBER          
C     COUT     - OUTPUT BLOCK
C     KNEXT    - POINTER TO NEXT LOCATION
C     IENT     - SWITCH  IENT=1  NORMAL
C                        IENT=2  END; OUTPUT LAST BLOCK
C                .      .    .                                       .
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     IREPTS   - NUMBER OF REPORTS             
C
C REMARKS: LIST CAVEATS, OTHER HELPFUL HINTS OR INFORMATION
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY4
C
C$$$
C
      SUBROUTINE DOSBYT(IARR, KRET, IUNO, COUT, KNEXT, IREPTS, IENT)
C
      CHARACTER*8   COUT (512)
      CHARACTER*8   CEND, CSTR, CBLK, CENDFL  
C
      INTEGER*8     IARR (300)
C
      DATA CEND    /'END_REPO'/
      DATA CSTR    /'STR_REPO'/
      DATA CBLK    /'        '/
      DATA CENDFL  /'ENDOFILE'/
      DATA KLIM    /512/
      DATA SHFTL   /4096/
C
C     THIS SUBROUTINE USES SBYTES TO FORM THE OUTPUT.
C     THE BLOCK IS WRITTEN WHEN FULL.
C  
C     IF IENT = 1   -   LOAD REPORT INTO BLOCK
C        IENT = 2   -   END CONDITION. WRITE LAST BLOCK.
C
      IF (IENT .EQ. 1) THEN
         IREPTS = IREPTS + 1 
         INEXT = KNEXT + KRET/2
         COUT(KNEXT) = CSTR
         IARR(3) = (INEXT*2 + 1)*SHFTL
         KNEXT = KNEXT + 1
C        PRINT 100,IARR(11),IARR(12),INEXT
 100     FORMAT(' STATION ',a8,a8,' INEXT IS ',i8)
         IF (INEXT .GE. KLIM) THEN
            INEXT = INEXT - KLIM
            IARR(3) = (INEXT*2 + 1)*SHFTL
            IST = 2*(KLIM - KNEXT + 1)
            CALL SBYTES(COUT(KNEXT), IARR, 0, 32, 0, IST)
            WRITE(IUNO) COUT
            KRET = KRET - IST 
            DO KK = 1,KLIM
               COUT(KK) = CBLK
            END DO
            KNEXT = 1
            IF (KRET .NE. 0) THEN
               CALL SBYTES(COUT, IARR(IST+1), 0, 32, 0, KRET)
               KNEXT = KRET/2 + 1
            END IF
            COUT(KNEXT) = CEND
            KNEXT = KNEXT + 1
          ELSE 
            CALL SBYTES(COUT(KNEXT), IARR, 0, 32, 0, KRET)
            KNEXT = KNEXT + KRET/2
            IF (KNEXT .GT. KLIM) THEN
  		WRITE(IUNO) COUT
		KNEXT = 1
		DO KK = 1,KLIM
		   COUT(KK) = CBLK
		END DO
	    END IF
            COUT(KNEXT) = CEND
            KNEXT = KNEXT + 1
            IF (KNEXT .GT. KLIM) THEN
               WRITE(IUNO) COUT
               KNEXT = 1
               DO KK = 1,KLIM
                  COUT(KK) = CBLK
               END DO
            END IF
          END IF
       ELSE
C
C         OUTPUT THE LAST BLOCK OF DATA
C
          IF (KNEXT.NE.1) THEN
             COUT(KNEXT) = CENDFL
             WRITE (IUNO) COUT
          ELSE
C
C            OUTPUT AN END-OF-FILE RECORD, FOR THE CASE
C            WHEN THERE IS NO DATA IN THE LAST BLOCK EXCEPT
C            THE END OF FILE INDICATOR
C
             COUT(1) = CENDFL
             DO  K=1,KLIM
                COUT(K) = CBLK
             END DO
             WRITE(IUNO) COUT
          END IF
       END IF
       RETURN
C
       END
