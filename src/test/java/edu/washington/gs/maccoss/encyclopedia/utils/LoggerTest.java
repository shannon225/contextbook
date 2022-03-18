package edu.washington.gs.maccoss.encyclopedia.utils;

import junit.framework.TestCase;

public class LoggerTest extends TestCase {
    Throwable cause() {
        return new RuntimeException("Cause.");
    }

    Throwable exception() {
        return new RuntimeException("Exception.", cause());
    }

    /**
     * No assertions; check the console output!
     * We expect to see both the exception stacktrace AND it's cause.
     */
    public void testLogException() {
        Logger.logException(exception());
    }

    public void testErrorException() {
        Logger.errorException(exception());
    }
}