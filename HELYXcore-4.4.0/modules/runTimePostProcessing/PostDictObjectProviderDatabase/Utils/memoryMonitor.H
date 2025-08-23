/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */
#pragma once

#include <cstddef>

class MemoryMonitor
{
public:
    /**
     * Returns the peak (maximum so far) resident set size (physical
     * memory use) measured in bytes, or zero if the value cannot be
     * determined on this OS.
     */
    static size_t getPeakRSS();
    static size_t getPeakRSSMByte() {return getPeakRSS() / (1024 * 1024);};

    /**
     * Returns the current resident set size (physical memory use) measured
     * in bytes, or zero if the value cannot be determined on this OS.
     */
    static size_t getCurrentRSS();
    static size_t getCurrentRSSMByte() {return getCurrentRSS() / (1024 * 1024);};
};
