#!/usr/bin/env bash
set -euo pipefail

#---------------------------------------------
# run folder 위치에서 수행
# cd /ruby/Analysis/Project/SMRT2_data/reDemulti_RUN/revio_01/r84135_20251003_081955/
# bash /ruby/Analysis/Project/SMRT2_data/reDemulti_RUN/revio_01/r84135_20251003_081955/test.py
#---------------------------------------------
# 설정
#---------------------------------------------
# 스레드 수: SGE(NSLOTS) > 환경변수 T > 16
T="${NSLOTS:-${T:-16}}"

# samtools 경로(환경에 맞게 조정). 기본값은 PATH상의 samtools
SAMTOOLS="${SAMTOOLS:-/ruby/Users/yeongjae0420/00.Tools/SAMTOOLS/samtools-1.19.2/samtools}"
if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
  SAMTOOLS=samtools
fi

# 출력 폴더
OUT="out_by_well"
mkdir -p "$OUT"

# 안전 옵션: 글롭 미매치 시 루프 건너뛰기
shopt -s nullglob

echo "[INFO] T=$T  SAMTOOLS=$SAMTOOLS"
echo "[INFO] OUT=$OUT"

#---------------------------------------------
# 1단계) 웰별 readstats 누적 생성
#   - 여러 BAM(동일 웰) → 한 개의 웰별 readstats로 append
#   - 필터: -F 0x900 (secondary/supplementary 제외)
#   - 태그: rq 필수, np 없으면 0으로 보정
#   - 출력: <length>\t<rq>\t<np>
#---------------------------------------------
# 기존 산출물 초기화(있다면 삭제)
rm -f "${OUT}"/*.readstats

for BAM in */hifi_reads/*.hifi_reads.bam; do
  WELL="$(echo "$BAM" | cut -d/ -f1)"
  PREFIX="${OUT}/${WELL}"

  echo "[INFO] Processing WELL=${WELL}"
  echo "       BAM=${BAM}"
  echo "       OUT=${PREFIX}.*"

  # append 모드로 누적
  "$SAMTOOLS" view -@ "$T" -F 0x900 "$BAM" \
  | awk 'BEGIN{OFS="\t"}{
      len=length($10); rq=""; np="";
      for(i=12;i<=NF;i++){
        split($i,t,":");
        if(t[1]=="rq") rq=t[3];
        else if(t[1]=="np") np=t[3];
      }
      # rq가 있을 때만 출력, np 없으면 0으로 보정
      if(rq!=""){
        if(np=="") np=0;
        print len, rq, np;
      }
    }' >> "${PREFIX}.readstats"
done

#---------------------------------------------
# 2단계) 파생 파일 생성
#   - _HiFi_Length.txt
#   - .hifi_reads_Phred.txt (per-read Q = -10*log10(1-acc), acc=1.0은 1-ε로 클램프)
#---------------------------------------------
for RS in "${OUT}"/*.readstats; do
  [ -s "$RS" ] || continue
  base="${RS%.readstats}"

  # 길이만
  cut -f1 "$RS" > "${base}_HiFi_Length.txt"

  # per-read Phred
  awk '{
    acc=$2; eps=1e-12; if(acc>=1.0){acc=1.0-eps}
    q=-10*log(1-acc)/log(10);
    printf("%.4f\n", q);
  }' "$RS" > "${base}.hifi_reads_Phred.txt"
done

#---------------------------------------------
# 3단계) 요약 통계(_stat.txt) - gawk 히스토그램 방식(메모리 절약)
#   출력 헤더:
#   Polymerase_Read_Bases, Polymerase_Reads, Polymerase_N50, Ave_Polymerase_Length,
#   Subreads_Bases, Subreads_Reads, Ave_Subreads_Length, Subreads_N50,
#   HiFi_Bases, HiFi_Reads, Ave_HiFi_Length, HiFi_N50, HiFi_MeanPass, HiFi_Quality
#   - HiFi_Quality: per-read Q의 "중앙값(median)"을 정수 Q로 출력
#   - HiFi_MeanPass: np의 산술평균을 반올림하여 정수로 출력
#---------------------------------------------
for RS in "${OUT}"/*.readstats; do
  base="${RS%.readstats}"
  STAT="${base}_stat.txt"

  if [ ! -s "$RS" ]; then
    # 빈 파일: 0/헤더만 출력
    {
      echo -e "Polymerase_Read_Bases\tPolymerase_Reads\tPolymerase_N50\tAve_Polymerase_Length\tSubreads_Bases\tSubreads_Reads\tAve_Subreads_Length\tSubreads_N50\tHiFi_Bases\tHiFi_Reads\tAve_HiFi_Length\tHiFi_N50\tHiFi_MeanPass\tHiFi_Quality"
      echo -e "0\t0\t0\t0\t-\t-\t-\t-\t0\t0\t0\t0\t0\tQ0"
    } > "$STAT"
    continue
  fi

  gawk 'BEGIN{OFS="\t"}
  {
    L=$1+0; A=$2; P=$3;           # A와 P는 결측 가능
    # 길이 통계(HiFi_Bases/Reads/N50/평균길이) - readstats가 rq 있는 리드만 포함하는 전제에서 계산
    cL[L] += 1;         # 길이별 개수(히스토그램)
    sumL  += L;         # HiFi_Bases
    nL    += 1;         # HiFi_Reads

    # passes 평균(존재하는 경우만 집계)
    if(A != "" ){              # 품질 히스토그램용(Q median)
      acc = A + 0.0;
      eps = 1e-12; if(acc>=1.0){acc=1.0-eps}
      Qf = -10*log(1-acc)/log(10);
      Qi = int(Qf + 0.5);      # 정수 Q로 반올림
      if(Qi<0) Qi=0;
      cQ[Qi] += 1;
      nQ     += 1;
    }
    if(P != ""){               # np 평균(존재하는 경우만)
      sumP += (P+0);
      nP   += 1;
    }
  }
  END{
    # 헤더
    print "Polymerase_Read_Bases","Polymerase_Reads","Polymerase_N50","Ave_Polymerase_Length",
          "Subreads_Bases","Subreads_Reads","Ave_Subreads_Length","Subreads_N50",
          "HiFi_Bases","HiFi_Reads","Ave_HiFi_Length","HiFi_N50","HiFi_MeanPass","HiFi_Quality" > "'"$STAT"'";

    if(nL==0){
      print 0,0,0,0,"-","-","-","-",0,0,0,0,0,"Q0" >> "'"$STAT"'";
      exit
    }

    # N50 (길이 히스토그램 내림차순 누적)
    half = sumL/2.0;
    asorti(cL, idxL, "@ind_num_desc");
    run=0; N50=0;
    for(i=1;i<=length(idxL);i++){
      L = idxL[i];
      run += L * cL[L];
      if(run >= half){ N50=L; break }
    }

    aveLen = int(sumL / nL);

    # Mean Pass = 산술평균 반올림(존재하는 np에 대해)
    meanP = 0;
    if(nP>0){
      meanP = int( (sumP / nP) + 0.5 );
    }

    # Q(median) = per-read 정수 Q의 중앙값
    Qmed = 0;
    if(nQ>0){
      mid = int((nQ+1)/2);          # 중앙 인덱스(하위 중앙)
      asorti(cQ, idxQ, "@ind_num_asc");  # Q 오름차순
      cum=0;
      for(i=1;i<=length(idxQ);i++){
        qv = idxQ[i];
        cum += cQ[qv];
        if(cum >= mid){ Qmed = qv; break }
      }
    }

    print 0,0,0,0,"-","-","-","-",sumL,nL,aveLen,N50,meanP,"Q"Qmed >> "'"$STAT"'";
  }' "$RS"
done

echo "[DONE] Per-well files are under: ${OUT}/"
