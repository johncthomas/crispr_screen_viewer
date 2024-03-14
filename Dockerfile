# syntax=docker/dockerfile:1

# Comments are provided throughout this file to help you get started.
# If you need more help, visit the Dockerfile reference guide at
# https://docs.docker.com/engine/reference/builder/
ARG PYTHON_VERSION=3.12
ARG DATASOURCE=./src/crispr_screen_viewer/data/test_db

FROM python:${PYTHON_VERSION}-slim as base

EXPOSE 8050

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# install git, gcc and g++
RUN --mount=target=/var/lib/apt/lists,type=cache,sharing=locked \
    --mount=target=/var/cache/apt,type=cache,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean
RUN apt-get update
RUN apt-get -y install gcc g++
RUN apt-get -y install git
RUN rm -rf /var/lib/apt/lists/*

# install pip requirements
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

RUN python -m pip install gunicorn

ARG DIR=/app
WORKDIR /app
ENV PATH="$PATH:/app/src/crispr_screen_viewer"

COPY . .
COPY ${DATASOURCE}/* /app/src/crispr_screen_viewer/data/test_db

RUN python -m pip install ${DIR}
