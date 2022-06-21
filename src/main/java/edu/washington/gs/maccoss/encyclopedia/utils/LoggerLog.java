package edu.washington.gs.maccoss.encyclopedia.utils;

import org.apache.commons.logging.Log;

import java.util.Objects;

public class LoggerLog implements Log {
	@Override
	public boolean isTraceEnabled() {
		return false;
	}

	@Override
	public boolean isDebugEnabled() {
		return false;
	}

	@Override
	public boolean isInfoEnabled() {
		return true;
	}

	@Override
	public boolean isWarnEnabled() {
		return true;
	}

	@Override
	public boolean isErrorEnabled() {
		return true;
	}

	@Override
	public boolean isFatalEnabled() {
		return true;
	}

	@Override
	public void trace(Object o) {
		trace(o, null);
	}

	@Override
	public void trace(Object o, Throwable throwable) {
		if (isTraceEnabled()) {
			Logger.logLine(Objects.toString(o));
			if (null != throwable) {
				Logger.logException(throwable);
			}
		}
	}

	@Override
	public void debug(Object o) {
		debug(o, null);
	}

	@Override
	public void debug(Object o, Throwable throwable) {
		if (isDebugEnabled()) {
			Logger.logLine(Objects.toString(o));
			if (null != throwable) {
				Logger.logException(throwable);
			}
		}
	}

	@Override
	public void info(Object o) {
		info(o, null);
	}

	@Override
	public void info(Object o, Throwable throwable) {
		if (isInfoEnabled()) {
			Logger.logLine(Objects.toString(o));
			if (null != throwable) {
				Logger.logException(throwable);
			}
		}
	}

	@Override
	public void warn(Object o) {
		warn(o, null);
	}

	@Override
	public void warn(Object o, Throwable throwable) {
		if (isWarnEnabled()) {
			Logger.errorLine(Objects.toString(o));
			if (null != throwable) {
				Logger.errorException(throwable);
			}
		}
	}

	@Override
	public void error(Object o) {
		error(o, null);
	}

	@Override
	public void error(Object o, Throwable throwable) {
		if (isErrorEnabled()) {
			Logger.errorLine(Objects.toString(o));
			if (null != throwable) {
				Logger.errorException(throwable);
			}
		}
	}

	@Override
	public void fatal(Object o) {
		fatal(o, null);
	}

	@Override
	public void fatal(Object o, Throwable throwable) {
		if (isFatalEnabled()) {
			Logger.errorLine(Objects.toString(o));
			if (null != throwable) {
				Logger.errorException(throwable);
			}
		}
	}
}
