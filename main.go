package main

import (
	"bytes"
	"fmt"
	"github.com/mjibson/go-dsp/dsputils"
	"github.com/mjibson/go-dsp/spectral"
	"github.com/mjibson/go-dsp/wav"
	"io"
	"log"
	"math"
	"os"
)

const winDenom = 2

type hit struct {
	t    float64
	freq float64
	pow  float64
}

type group struct {
	freqs []float64
}

func main() {
	if err := stuff(); err != nil {
		log.Printf("Error: %v\n", err)
	}
}

func stuff() error {
	hits, err := getHits()
	if err != nil {
		return err
	}

	hits = mergeHits(hits)
	for _, h := range hits {
		fmt.Printf("%v: %f Hz (%f)\n", fmtSeconds(h.t), h.freq, h.pow)
	}

	groups := findGroups(hits)
	for _, g := range groups {
		fmt.Printf("%v\n", g.String())
	}

	return nil
}

func getHits() ([]hit, error) {
	w, err := wav.New(os.Stdin)
	if err != nil {
		return nil, err
	}

	log.Printf("format: %d, channels: %d, sample rate: %d, byte rate: %d, bps: %d, samples: %d, duration: %v\n",
		w.Header.AudioFormat, w.Header.NumChannels, w.Header.SampleRate,
		w.Header.ByteRate, w.Header.BitsPerSample, w.Samples, w.Duration)

	winSize := int(w.Header.SampleRate * uint32(w.Header.NumChannels) / winDenom)

	cursor := 0
	hits := make([]hit, 0)

	for {
		rawSamples, err := w.ReadFloats(winSize)

		if err != nil {
			if err == io.EOF || err == io.ErrUnexpectedEOF {
				break
			} else {
				return nil, err
			}
		}

		if len(rawSamples) < winSize {
			break
		}

		samples := flattenChannels(int(w.Header.NumChannels), rawSamples)
		pxx, freqs := spectral.Pwelch(samples, float64(w.Header.SampleRate), &spectral.PwelchOptions{NFFT: 16384, Scale_off: true})
		maxPower, maxPowerFreq := findPeak(pxx, freqs, 32, 880)
		pMaxPower, pMaxPowerFreq := findPeak(pxx, freqs, 55, 170)

		harmonic := true

		if !dsputils.Float64Equal(maxPowerFreq, pMaxPowerFreq) {
			if maxPowerFreq < pMaxPowerFreq || !areHarmonic(maxPower, pMaxPower) {
				fmt.Printf("t = %v: %f !~ %f\n", fmtSeconds(float64(cursor)/float64(w.Header.SampleRate)), maxPowerFreq, pMaxPowerFreq)
				harmonic = false
			}
		}

		v1Freq := pMaxPowerFreq * 1.1
		v1Pow := powerAtFreq(v1Freq, pxx, freqs)

		p1Freq := pMaxPowerFreq * 1.2
		p1Pow := powerAtFreq(p1Freq, pxx, freqs)

		v2Freq := pMaxPowerFreq * 1.3
		v2Pow := powerAtFreq(v2Freq, pxx, freqs)

		p2Freq := pMaxPowerFreq * 1.4
		p2Pow := powerAtFreq(p2Freq, pxx, freqs)

		if harmonic && pMaxPower >= 0.00 && pMaxPowerFreq >= 55.0 && pMaxPowerFreq <= 160.0 && v1Pow < p1Pow && v2Pow < p2Pow {
			hits = append(hits, hit{t: float64(cursor) / float64(w.Header.SampleRate), freq: pMaxPowerFreq, pow: pMaxPower})
		}

		cursor += len(rawSamples) / int(w.Header.NumChannels)
	}

	return hits, nil
}

func flattenChannels(numChannels int, samples []float32) []float64 {
	outLen := len(samples) / numChannels
	out := make([]float64, outLen)
	for i := 0; i < outLen; i++ {
		var agg float64 = 0
		for j := 0; j < numChannels; j++ {
			agg += float64(samples[(i*numChannels)+j])
		}
		out[i] = agg / float64(numChannels)
	}

	return out
}

func findPeak(pxx, freqs []float64, min, max float64) (float64, float64) {
	var maxPower float64 = 0
	var maxPowerFreq float64 = 0

	for i, p := range pxx {
		freq := freqs[i]

		if freq < min || freq > max {
			continue
		}

		if p > maxPower {
			maxPower = p
			maxPowerFreq = freq
		}
	}

	return maxPower, maxPowerFreq
}

func powerAtFreq(freq float64, pxx, freqs []float64) float64 {
	for i := 0; i+1 < len(freqs); i++ {
		if freq >= freqs[i] && freq <= freqs[i+1] {
			return pxx[i]
		}
	}

	return 0
}

func fmtSeconds(seconds float64) string {
	min := int(math.Floor(seconds / 60))
	sec := int(math.Floor(seconds)) % 60
	dsec := int(math.Floor(10 * (seconds - math.Floor(seconds))))

	return fmt.Sprintf("%02d:%02d.%d", min, sec, dsec)
}

func areHarmonic(a, b float64) bool {
	min := math.Min(a, b)
	max := math.Max(a, b)

	r := max / min
	i0 := (r / 0.4)
	i1 := math.Floor(i0 + 0.5)
	id := math.Abs(i1 - i0)
	return id <= 0.05
}

func mergeHits(src []hit) []hit {
	out := make([]hit, 0)

	for i := 0; i < len(src); i++ {
		base := src[i]
		last := base
		grp := make([]hit, 1)
		grp[0] = base

		for j := i + 1; j < len(src); j++ {
			top := src[j]
			if top.t-last.t <= 1.1 && math.Abs(top.freq-last.freq) < 4.0 {
				grp = append(grp, top)
				last = top
			} else {
				break
			}
		}

		if len(grp) == 1 {
			out = append(out, base)
		} else {
			weightedFreqAcc := 0.0
			powerAcc := 0.0
			maxPower := 0.0

			for _, h := range grp {
				weightedFreqAcc += h.freq * h.pow
				powerAcc += h.pow
				if h.pow > maxPower {
					maxPower = h.pow
				}
			}

			out = append(out, hit{t: base.t, freq: weightedFreqAcc / powerAcc, pow: maxPower})
			i += len(grp) - 1
		}
	}

	return out
}

func findGroups(hits []hit) []group {
	out := make([]group, 0)

	for i := 0; i < len(hits); i++ {
		base := hits[i]
		last := base
		grp := make([]hit, 1)
		grp[0] = base

		consumed := 0
		for j := i + 1; j < len(hits); j++ {
			top := hits[j]
			tDiff := top.t - last.t
			if tDiff < 2.9 {
				consumed++
				continue
			}

			if tDiff <= 6.1 {
				if !areFreqsGroupable(base.freq, top.freq) {
					consumed++
					continue
				}

				consumed++
				grp = append(grp, top)
				last = top
			} else {
				break
			}
		}

		freqs := make([]float64, len(grp))
		for idx, hit := range grp {
			freqs[idx] = hit.freq
		}

		sgrp := group{freqs: freqs}
		out = append(out, sgrp)

		i += consumed
	}

	return out
}

func areFreqsGroupable(f0, f1 float64) bool {
	return areFreqsSimilar(f0, f1) || areFreqsRelated(f0, f1)
}

func areFreqsSimilar(f0, f1 float64) bool {
	r := f0 / f1
	return r >= 0.9 && r <= 1.1
}

func areFreqsRelated(f0, f1 float64) bool {
	return (math.Abs((f0/f1)-1.5) <= 0.2) || (math.Abs((f1/f0)-1.5) <= 0.2)
}

func (g *group) String() string {
	b := &bytes.Buffer{}

	for idx, freq := range g.freqs {
		if idx > 0 {
			b.WriteString(", ")
		}

		b.WriteString(fmt.Sprintf("%.0f", freq))
	}

	return b.String()
}
